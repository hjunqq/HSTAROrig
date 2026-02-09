    subroutine modf_element_local_direction  !20221102

    integer(ink) igroup,index,nnode,ielgroup,ielem,inode,idimn
    integer(ink),pointer::lnods(:)

    real(irk) djacb
    real(irk),pointer::elcod(:,:)
    real(irk),allocatable::a3(:)

    do igroup=1,ngroup
        index = group(igroup)%index
        nnode = elkn(index)%el_field(1)%nnode_f

        if(index/=20.or.index/=21.or.index/=22) cycle
        allocate(a3(ndimn))
        DO ielgroup = 1,group(igroup)%nelgroup
            ielem = group(igroup)%list(ielgroup)

            lnods=>element(ielem)%field(1)%lnods_f
            elcod=>element(ielem)%field(1)%elcod_f

            if(index==22)then
                call normal_local_plate(index,ndimn,elcod,a3)
                call direct_p4(ndimn,a3,element(ielem)%rotation,element(ielem)%point_direct,coord)
            else
                djacb =sqrt(sum((elcod(1:ndimn,2)-elcod(1:ndimn,1))**2))
                a3(:)=(elcod(:,2)-elcod(:,1))/djacb
                call direct_beam(ndimn,a3,element(ielem)%rotation,element(ielem)%point_direct,coord)
            endif

            do inode=1,nnode
                do idimn=1,ndimn
                    elcod(idimn,inode)=element(ielem)%rotation(idimn,:).d.coord(:,lnods(inode))
                end do
            end do
            element(ielem)%field(1)%elcod_f=elcod
            nullify(elcod,lnods)
        end do
        deallocate(a3)
    end do
    end subroutine modf_element_local_direction !20221102

    subroutine dudx  !20220108
    integer(ink) ivar,igroup,matno,ielgroup,ielem,nevab,itotv,iieq,  &
        ivalue,iintf,nintf,inode,idofn,jnode,   &
        iblks_i,iincs_i,istep_i,ivalue_point,jvar,jtotv,njntf
    real   (irk) e
    real   (irk),allocatable::stfor_bar(:),eload(:),value(:)
    real   (irk),pointer::fstif(:,:),rintf(:)
    integer(ink),pointer::ldofs(:),listf(:)

    !do ivalue=1,mvalue
    !Value_observ(ivalue)%dudx=0.
    !end do
    allocate(stfor_bar(ntotv))
    do ivar=1,npara
        write(7,*)'ivar=',ivar,'para_back(ivar)%name=',para_back(ivar)%name
        stfor_bar=0.
        if(para_back(ivar)%name=='E'.or.para_back(ivar)%name=='PERM')then  !20230430

            DO igroup =1,ngroup  !igroup
                if (appear(igroup)<=0) cycle
                matno = group(igroup)%matno
                if(para_back(ivar)%name=='E')then  !20230430
                    if(props(matno)%mechanical%solid%ie/=ivar)cycle
                    e=xvalue(ivar)
                    if(para_back(ivar)%mode_transform==0)then
                        e=xvalue(ivar)/para_back(ivar)%factor
                    elseif(para_back(ivar)%mode_transform==1)then
                        e=1./xvalue(ivar)/para_back(ivar)%factor
                    endif
                    !write(7,*)'igroup=',igroup,'e=',e
                endif
                if(para_back(ivar)%name=='PERM')then !20230430
                    if(props(matno)%mechanical%fluid%iperm/=ivar)cycle
                    e=xvalue(ivar)
                    if(para_back(ivar)%mode_transform==0)then
                        e=xvalue(ivar)/para_back(ivar)%factor
                    elseif(para_back(ivar)%mode_transform==1)then
                        e=1./xvalue(ivar)/para_back(ivar)%factor
                    endif
                endif

                nevab=size(element(group(igroup)%list(1))%field(1)%ldofs_f)
                allocate(eload(nevab),value(nevab))
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    fstif=>element(ielem)%field(1)%khandmc(1)%fstif
                    ldofs=>element(ielem)%field(1)%ldofs_f
                    if (Bparameter==1) value=result_zero(ldofs)
                    if (Bparameter==2) value=deltafi(ldofs)
                    !fstif=fstif/e
                    eload=MATMUL(fstif/e,value)
                    stfor_bar(ldofs)=stfor_bar(ldofs)-eload
                    nullify(fstif,ldofs)
                end do
                deallocate(eload,value)
            end do  !igroup

            !!
            rvector=0.
            do itotv=1,ntotv
                if (totveq(itotv)==0)cycle
                rvector(totveq(itotv))=rvector(totveq(itotv))+stfor_bar(itotv)
            end do
            do itotv=1,ntotv
                nintf=trans(itotv)%nintf
                if (nintf==0)cycle
                do iintf=1,nintf
                    iieq=totveq(trans(itotv)%listf(iintf))
                    if(iieq/=0)rvector(iieq)=rvector(iieq)+stfor_bar(itotv)*trans(itotv)%rintf(iintf)
                end do
            end do
            !!
            operation='SOLVE'
            call solve
            !! 找出观测点处的du/dx

            do ivalue=1,mvalue
                if(Value_observ(ivalue)%ic==0)cycle
                iblks_i=Value_observ(ivalue)%iblks
                iincs_i=Value_observ(ivalue)%iincs
                istep_i=Value_observ(ivalue)%istep
                idofn =lmdofn(Value_observ(ivalue)%idofn)

                ivalue_point=Value_observ(ivalue)%ivalue_point
                if(iblks_i==iblks.and.iincs_i==iincs.and.istep_i==istep)then
                    nintf=para_points(ivalue_point)%nintf
                    listf=>para_points(ivalue_point)%listf
                    rintf=>para_points(ivalue_point)%rintf
                    Value_observ(ivalue)%dudx(ivar)=dot_product(rintf,result(nodfn(idofn,listf)))
                    nullify(listf,rintf)
                endif  ! if(iblks_i==iblks....)
            end do  !ivalue
            !write(7,*)'istep=',istep,'ivar=',ivar,'value_observe%dudx='
            !write(7,*)Value_observ(:)%dudx(ivar)
        endif  !20230430
        !!!!
    end do  !ivar

    deallocate(stfor_bar)
    end subroutine dudx !20220108

    subroutine cmatrix_c_formation  !20210328
    integer(ink) igapb,igroup,nline_g_sc,npairs_sc,ipoin,jpoin,jdimn,jtotv,nintf,  &
        iintf,il0,kdimn,kpoin,iieq
    real   (irk) dislocal_s
    real   (irk),allocatable::unitg(:)
    real   (irk),pointer::Ks(:,:),Kc(:,:)

    allocate(unitg(ndimn))

    !write(7,*)'cmatrix_c='
    do igapb=1,nrcsteel
        igroup=rc_steel(igapb)%listgroup_c
        if(appear_process(igroup,iblks)==0)cycle
        nline_g_sc=rc_steel(igapb)%nline_g_sc
        do il0=1,nline_g_sc
            !write(7,*)'steel_group=',igapb,'line_order=',il0
            npairs_sc=rc_steel(igapb)%line_g_sc(il0)%npairs_sc

            rc_steel(igapb)%line_g_sc(il0)%cmatrix_c=0.

            do ipoin=1,npairs_sc
                rvector=0.
                jpoin=rc_steel(igapb)%line_g_sc(il0)%pairnode_sc(ipoin)

                do jdimn=1,ndimn
                    jtotv=nodfn(jdimn,jpoin)
                    nintf=trans(jtotv)%nintf
                    if(nintf/=0) then
                        iieq=totveq(jtotv)
                        if(iieq/=0)rvector(iieq)=0.
                        do iintf=1,nintf
                            iieq=totveq(trans(jtotv)%listf(iintf))
                            if(iieq/=0) &
                                rvector(iieq)=rvector(iieq)+rc_steel(igapb)%line_g_sc(il0)%rot_sc(jdimn,ipoin)*trans(jtotv)%rintf(iintf)
                        end do
                    else
                        if(totveq(jtotv)/=0) &
                            rvector(totveq(jtotv))=rvector(totveq(jtotv))+rc_steel(igapb)%line_g_sc(il0)%rot_sc(jdimn,ipoin)

                    endif
                end do !jdimn

                operation='SOLVE'
                call solve

                do kpoin=1,npairs_sc
                    jpoin=rc_steel(igapb)%line_g_sc(il0)%pairnode_sc(kpoin)
                    unitg=0.
                    do jdimn=1,ndimn
                        jtotv=nodfn(jdimn,jpoin)
                        nintf=trans(jtotv)%nintf
                        if(nintf/=0) then
                            do iintf=1,nintf
                                itotv=trans(jtotv)%listf(iintf)
                                unitg(jdimn)=unitg(jdimn)+result(itotv)*trans(jtotv)%rintf(iintf)
                            end do
                        else
                            unitg(jdimn)=result(jtotv)
                        endif
                    end do !jdimn
                    dislocal_s=dot_product(unitg,rc_steel(igapb)%line_g_sc(il0)%rot_sc(:,kpoin))
                    rc_steel(igapb)%line_g_sc(il0)%cmatrix_c(kpoin,ipoin)=dislocal_s

                end do  !kpoin
            end do  !ipoin


            Ks=>rc_steel(igapb)%line_g_sc(il0)%kmatrix_s
            Kc=>rc_steel(igapb)%line_g_sc(il0)%cmatrix_c
            rc_steel(igapb)%line_g_sc(il0)%ikscr=(Ks.x.Kc)
            do ipoin=1,npairs_sc
                rc_steel(igapb)%line_g_sc(il0)%ikscr(ipoin,ipoin)=   &
                    rc_steel(igapb)%line_g_sc(il0)%ikscr(ipoin,ipoin)+1.
            end do

            !write(7,*)'Ks='
            !       do ipoin=1,npairs_sc
            !   write(7,10)Ks(ipoin,:)
            !       end do
            !
            !   write(7,*)'Kc='
            !       do ipoin=1,npairs_sc
            !   write(7,10)Kc(ipoin,:)
            !       end do
            !            write(7,*)'ikscr='
            !    do ipoin=1,npairs_sc
            !write(7,10)rc_steel(igapb)%line_g_sc(il0)%ikscr(ipoin,:)
            !    end do
            nullify(ks,kc)

        end do !il0

    end do  !igapb
    deallocate(unitg)
10  format(15e15.5)

    end subroutine cmatrix_c_formation  !20210328

    subroutine Tcmatrix_c_formation  !20210411
    integer(ink) igapb,igroup,nline_g_w,npairs_wc,ipoin,jpoin,jdimn,jtotv,nintf,  &
        iintf,il0,kdimn,kpoin,iieq
    real   (irk) dislocal_s,coef
    real   (irk),pointer::Ks(:,:),Kc(:,:)


    !write(7,*)'cmatrix_c='
    coef=theta1*ditime
    do igapb=1,nwcpipe
        igroup=wc_pipe(igapb)%listgroup_c
        if(appear_process(igroup,iblks)==0)cycle
        nline_g_w=wc_pipe(igapb)%nline_g_w

        do il0=1,nline_g_w
            !write(7,*)'steel_group=',igapb,'line_order=',il0
            npairs_wc=wc_pipe(igapb)%line_g_w(il0)%npairs_wc

            wc_pipe(igapb)%line_g_w(il0)%cmatrix_c=0.

            do ipoin=1,npairs_wc
                rvector=0.
                jpoin=wc_pipe(igapb)%line_g_w(il0)%pairnode_wc(ipoin)
                jtotv=nodfn(lmdofn(10),jpoin)
                nintf=trans(jtotv)%nintf
                if(nintf/=0) then
                    iieq=totveq(jtotv)
                    if(iieq/=0)rvector(iieq)=0.
                    do iintf=1,nintf
                        iieq=totveq(trans(jtotv)%listf(iintf))
                        if(iieq/=0) &
                            rvector(iieq)=rvector(iieq)+1.*trans(jtotv)%rintf(iintf)
                    end do
                else
                    if(totveq(jtotv)/=0) &
                        rvector(totveq(jtotv))=rvector(totveq(jtotv))+1.

                endif


                operation='SOLVE'
                call solve

                do kpoin=1,npairs_wc
                    jpoin=wc_pipe(igapb)%line_g_w(il0)%pairnode_wc(kpoin)


                    jtotv=nodfn(lmdofn(10),jpoin)
                    nintf=trans(jtotv)%nintf
                    dislocal_s=0.
                    if(nintf/=0) then
                        do iintf=1,nintf
                            itotv=trans(jtotv)%listf(iintf)
                            dislocal_s=dislocal_s+result(itotv)*trans(jtotv)%rintf(iintf)
                        end do
                    else
                        dislocal_s=result(jtotv)
                    endif
                    wc_pipe(igapb)%line_g_w(il0)%cmatrix_c(kpoin,ipoin)=dislocal_s

                end do  !kpoin
            end do  !ipoin


            Ks=>wc_pipe(igapb)%line_g_w(il0)%kmatrix_w
            Kc=>wc_pipe(igapb)%line_g_w(il0)%cmatrix_c
            wc_pipe(igapb)%line_g_w(il0)%idcr=(Kc.x.Ks)
            do ipoin=1,npairs_wc
                wc_pipe(igapb)%line_g_w(il0)%idcr(ipoin,ipoin)=   &
                    coef*wc_pipe(igapb)%line_g_w(il0)%idcr(ipoin,ipoin)+1.
            end do

            !write(7,*)'Ks='
            !       do ipoin=1,npairs_wc
            !   write(7,10)Ks(ipoin,:)
            !       end do
            !
            !   write(7,*)'Kc='
            !       do ipoin=1,npairs_wc
            !   write(7,10)Kc(ipoin,:)
            !       end do
            !               write(7,*)'idcr='
            !       do ipoin=1,npairs_wc
            !   write(7,10)wc_pipe(igapb)%line_g_w(il0)%idcr(ipoin,:)
            !       end do
            nullify(ks,kc)

        end do !il0

    end do  !igapb

10  format(15e15.5)

    end subroutine Tcmatrix_c_formation  !20210411

    subroutine cmatrix_dtv_formation(cmatrix_dtv,inv_cmatrix_dtv2)  !20230216
    character(10)fieldid
    integer(ink) ifixset,idofix,jfixset,idofn,igroup,nrfields,ikh,index,nevab,ifield,ipoin,inode,nintf,mfixset
    real   (irk) coef,cmatrix_dtv(:,:),inv_cmatrix_dtv2(:,:)
    real   (irk),allocatable::value(:),qtemp(:),dtemp(:),interpt(:,:),observt(:,:)
    integer(ink),pointer::ldofs(:),listf(:),mlist(:)
    real   (irk),pointer::fstif(:,:),rintf(:)


    allocate(interpt(nfixsets,nfixsets),observt(nfixsets,nfixsets))
    cmatrix_dtv=0.


    allocate(dtemp(ntotv),qtemp(ntotv))

    !write(7,*)'cmatrix_c='

    do ifixset=1,nfixsets
        dtemp=0.
        qtemp=0.
        do idofix=1,ndofix
            idofn=prescrib(idofix)%ldofix
            mfixset=prescrib(idofix)%mfixset !20231130
            mlist=>prescrib(idofix)%mlist !20231130
            do i0=1,mfixset  !20231130
                jfixset=prescrib(idofix)%mlist(i0) !20231130
                if(jfixset/=ifixset)cycle
                dtemp(idofn)=1.*prescrib(idofix)%rintf(i0) !20231130
            end do  !20231130
            nullify(mlist) !20231130
        end do

        DO igroup =1,ngroup

            if (appear(igroup)>0) then
                ! get information from the group level

                nrfields=group(igroup)%nrfields
                fieldid=group(igroup)%fieldid
                if (nrfields/=1)     cycle
                if (fieldid/='T'.and.fieldid/='W')     cycle

                index  = group(igroup)%index
                nevab    =elkn(index)%el_field(1)%nnode_f
                allocate(value(nevab))

                do ifield=1,nrfields
                    do ikh=1,2

                        coef=1.
                        if (type_problem=='Q'.and.ikh==2) cycle
                        if (type_problem=='S'.and.ikh==1) coef=ditime    !theta1*ditime

                        DO ielgroup = 1,group(igroup)%nelgroup
                            ielem = group(igroup)%list(ielgroup)

                            if (associated(element(ielem)%field(ifield)%khandmc(ikh)%fstif)) then

                                fstif=>element(ielem)%field(ifield)%khandmc(ikh)%fstif
                                ldofs=>element(ielem)%field(ifield)%ldofs_f

                                value=coef*dtemp(ldofs)
                                ic=size(fstif,dim=2)
                                if (ic==1)then
                                    do idofn=1,nevab
                                        qtemp(ldofs(idofn))=qtemp(ldofs(idofn))+fstif(idofn,1)*value(idofn)
                                    end do
                                else
                                    qtemp(ldofs)=qtemp(ldofs)+MATMUL(fstif,value)
                                endif
                                nullify(fstif,ldofs)
                            endif

                        end do       !!ielgroup
                    end do        !!end do ikh
                end do     !! end do ifield
            endif
            deallocate(value)
        end do     !!  for igroup




        rvector=0.
        do itotv=1,ntotv
            if(totveq(itotv)/=0) &
                rvector(totveq(itotv))=rvector(totveq(itotv))-qtemp(itotv)
        end do


        operation='SOLVE'
        call solve
        !!!!!!!
        do ipoin=1, Npoints_pbx
            nintf=para_points(ipoin)%nintf
            listf=>para_points(ipoin)%listf
            rintf=>para_points(ipoin)%rintf
            cmatrix_dtv(ipoin,ifixset)=dot_product(rintf,result(listf))
            nullify(listf,rintf)
        end do



    end do  !ifixset

    interpt=transpose(cmatrix_dtv).x.cmatrix_dtv
    observt=0.
    do ifixset=1,nfixsets
        observt(ifixset,ifixset)=1.
    end do
    call householder(interpt,observt,inv_cmatrix_dtv2)


    deallocate(qtemp,dtemp)
    deallocate(interpt,observt)

    end subroutine cmatrix_dtv_formation  !20230216

    subroutine forAdirect_back_analysis !20150925

    integer(ink) igapb,npgblock,onetwo,ipoin,jpoin,ipair,idimn,itotvbt,itotv,ielem, &
        kpoin,lpoin,jdimn,jtotv,jtotvbt,nevab,ieqx,ievab,igaps,nnode,ii,matno,index,order_int,ngaus, &
        ielgroup,igaus,igroup,jgroup,igapbf,kdimn,icg,ipoin0
    real   (irk) coef,tvol,density,djacb,rr,fact,xij
    real   (irk),allocatable::eldis(:),eload(:),rot(:,:),disgi(:,:),disli(:,:),dist(:),center(:),mass_inertia(:),gpcod(:),uirigid(:,:),disgi0(:,:)
    real   (irk),allocatable::coordx(:)
    integer(ink),pointer::ldofs(:)
    real   (irk),pointer::estif(:,:)

    kdimn=ndimn
    if(block_stab==1)kdimn=3*(ndimn-1)
    allocate(center(ndimn),rot(ndimn,ndimn),dist(ndimn),disgi(ndimn,3*(ndimn-1)))
    if(block_stab/=1) allocate(disli(ndimn,3*(ndimn-1)))
    if(block_stab==1) allocate(disli(3*(ndimn-1),3*(ndimn-1)))
    if(type_problem=='F')allocate(mass_inertia(3*(ndimn-1)))

    write(7,*)'foradirect'

    allocate(gpcod(ndimn))
    do igapbf=1,nbackf
        igapb=backf(igapbf)%groupb
        if(gapb(igapb)%nrdof==0) cycle   !tcl
        gapb(igapb)%center=0.;
        if(type_problem=='F')gapb(igapb)%mass_inertia=0.
        tvol=0.
        do igroup=1,gapb(igapb)%ngroupb
            jgroup=gapb(igapb)%listgroupb(igroup)
            matno = group(jgroup)%matno
            index = group(jgroup)%index
            density=props(matno)%mechanical%solid%density
            order_int=elkn(index)%el_field(1)%order_intrules(1)
            ngaus = elkn(index)%ggaus(order_int)%ngaus
            do ielgroup=1,group(jgroup)%nelgroup
                ielem = group(jgroup)%list(ielgroup)

                do igaus=1,ngaus
                    djacb=element(ielem)%egaus(order_int)%djacb(igaus)
                    gpcod=element(ielem)%egaus(order_int)%gpcod(:,igaus)
                    gapb(igapb)%center=gapb(igapb)%center+djacb*gpcod*density
                    tvol=tvol+djacb*density
                end do
            enddo
        end do
        gapb(igapb)%center=gapb(igapb)%center/tvol
        write(7,*)'tvol=',tvol,'igapb=',igapb,'center=', gapb(igapb)%center

        if(type_problem=='F')then
            do igroup=1,gapb(igapb)%ngroupb
                jgroup=gapb(igapb)%listgroupb(igroup)
                matno = group(jgroup)%matno
                index = group(jgroup)%index
                density=props(matno)%mechanical%solid%density
                order_int=elkn(index)%el_field(1)%order_intrules(1)
                ngaus = elkn(index)%ggaus(order_int)%ngaus
                do ielgroup=1,group(jgroup)%nelgroup
                    ielem = group(jgroup)%list(ielgroup)

                    do igaus=1,ngaus
                        djacb=element(ielem)%egaus(order_int)%djacb(igaus)
                        gpcod=element(ielem)%egaus(order_int)%gpcod(:,igaus)
                        if(ndimn==2)then
                            rr=(gpcod(1)-gapb(igapb)%center(1))**2+(gpcod(2)-gapb(igapb)%center(2))**2
                            gapb(igapb)%mass_inertia(1)=gapb(igapb)%mass_inertia(1)+djacb*density
                            gapb(igapb)%mass_inertia(2)=gapb(igapb)%mass_inertia(2)+djacb*density
                            gapb(igapb)%mass_inertia(3)=gapb(igapb)%mass_inertia(3)+djacb*rr*density
                        elseif(ndimn==3)then
                            gapb(igapb)%mass_inertia(1)=gapb(igapb)%mass_inertia(1)+djacb*density
                            gapb(igapb)%mass_inertia(2)=gapb(igapb)%mass_inertia(2)+djacb*density
                            gapb(igapb)%mass_inertia(3)=gapb(igapb)%mass_inertia(3)+djacb*density
                            rr=(gpcod(2)-gapb(igapb)%center(2))**2+(gpcod(3)-gapb(igapb)%center(3))**2
                            gapb(igapb)%mass_inertia(4)=gapb(igapb)%mass_inertia(4)+djacb*rr*density
                            rr=(gpcod(1)-gapb(igapb)%center(1))**2+(gpcod(3)-gapb(igapb)%center(3))**2
                            gapb(igapb)%mass_inertia(5)=gapb(igapb)%mass_inertia(5)+djacb*rr*density
                            rr=(gpcod(2)-gapb(igapb)%center(2))**2+(gpcod(1)-gapb(igapb)%center(1))**2
                            gapb(igapb)%mass_inertia(6)=gapb(igapb)%mass_inertia(6)+djacb*rr*density
                        endif

                    end do
                enddo
            end do
            write(7,*)'igapb=',igapb,'gapb(igapb)%mass_inertia=', gapb(igapb)%mass_inertia
        endif
    enddo

    deallocate(gpcod)


    do igapbf=1,nbackf
        igapb=backf(igapbf)%groupb

        if(gapb(igapb)%nrdof==0) cycle   !tcl
        npgblock=gapb(igapb)%npgblock
        center=gapb(igapb)%center
        if(type_problem=='F')mass_inertia=gapb(igapb)%mass_inertia

        !write(7,*)'npdisp++++++++++++++'
        do jpoin=1,gapb(igapb)%npblock
            ipoin=gapb(igapb)%nodeblock(jpoin)
            dist=0.
            do jdimn=1,ndimn
                dist(jdimn)=coord(jdimn,ipoin)-center(jdimn)
            end do
            disgi=0.
            do idimn=1,ndimn
                disgi(idimn,idimn)=1.
            end do

            if(ndimn==2)then
                disgi(1,3)=-dist(2)
                disgi(2,3)= dist(1)
            elseif(ndimn==3)then
                disgi(1,5)= dist(3)
                disgi(1,6)=-dist(2)

                disgi(2,4)=-dist(3)
                disgi(2,6)= dist(1)

                disgi(3,4)= dist(2)
                disgi(3,5)=-dist(1)
            endif

            do idimn=1,gapb(igapb)%nrdof
                do jdimn=1,ndimn
                    gapb(igapb)%npdisp(jdimn,jpoin,idimn)=disgi(jdimn,idimn) !存每个点的位移，为计算A(T)F
                enddo
            end do
            if(block_stab==1)then
                do jdimn=ndimn+1,3*(ndimn-1)
                    gapb(igapb)%npdisp(jdimn,jpoin,jdimn)=1. !存每个点的位移，为计算A(T)F
                enddo
            endif
        enddo



        if(type_problem=='F')then
            do idimn=1,gapb(igapb)%nrdof
                gapb(igapb)%rstiff(idimn,idimn)=-mass_inertia(idimn) !20121216
            enddo !idimn
        endif

        if (restart_ctt/=0) cycle
        icg=0
        if(backf(igapbf)%mdism>npgblock*kdimn)then
            icg=1
            allocate(uirigid(backf(igapbf)%mdism,gapb(igapb)%nrdof))
        endif

        write(7,*)'icg=',icg


        do kpoin=1,backf(igapbf)%mdism
            allocate(coordx(ndimn))
            do idimn=1,ndimn
                coordx(idimn)=0.
                do i0=1,backf(igapbf)%relat(kpoin)%nintf
                    coordx(idimn)=coordx(idimn)+coord(idimn,backf(igapbf)%relat(kpoin)%listp(i0))*backf(igapbf)%relat(kpoin)%rintf(i0)
                end do
            end do
            !write(7,*)'kpoin=',kpoin,'coordx=',coordx
            dist=0.
            do jdimn=1,ndimn
                dist(jdimn)=coordx(jdimn)-center(jdimn)
            end do
            disgi=0.
            do idimn=1,ndimn
                disgi(idimn,idimn)=1.
            end do
            if(ndimn==2)then
                disgi(1,3)=-dist(2)
                disgi(2,3)= dist(1)
            elseif(ndimn==3)then
                disgi(1,5)= dist(3)
                disgi(1,6)=-dist(2)
                disgi(2,4)=-dist(3)
                disgi(2,6)= dist(1)
                disgi(3,4)= dist(2)
                disgi(3,5)=-dist(1)
            endif

            do idimn=1,gapb(igapb)%nrdof
                itotvbt=npgblock*kdimn+idimn
                jdimn=backf(igapbf)%listdim(kpoin)
                if(icg==0)gapb(igapb)%cmatrix(kpoin,itotvbt)=disgi(jdimn,idimn)
                if(icg==1)uirigid(kpoin,idimn)=disgi(jdimn,idimn)
            enddo !idimn
            deallocate(coordx)

        enddo !kpoin

        if(icg==1)then  !2015/11/28

            gapb(igapb)%cmatrix(1:npgblock*kdimn,npgblock*kdimn+1:npgblock*kdimn+gapb(igapb)%nrdof)=transpose(gapb(igapb)%uireact).x.uirigid

            deallocate(uirigid)
        endif   !2015/11/28


        do kpoin=1,npgblock
            onetwo=gapb(igapb)%nodegblock_onetwo(kpoin) !前四个点还是后四个点，也就是第一点还是第二点
            if(onetwo==1)coef=1.
            if(onetwo==2)coef=-1.
            igaps=gapb(igapb)%nodegblock_igaps(kpoin)
            ipair=gapb(igapb)%nodegblock_ipairs(kpoin)
            rot=gaps(igaps)%rot(:,:,ipair)


            dist=0.
            if(contactpe==1)then
                lpoin=gapb(igapb)%nodegblock(kpoin)
                do jdimn=1,ndimn
                    dist(jdimn)=coord(jdimn,lpoin)-center(jdimn)
                end do
            elseif(contactpe==2)then
                !nnode=2
                !if(ndimn==3)nnode=4
                nnode=size(gaps(igaps)%pairnode(:,ipair))/2  !2017/02/14

                do ii=1,nnode
                    jpoin=gaps(igaps)%pairnode((onetwo-1)*nnode+ii,ipair)
                    do jdimn=1,ndimn
                        dist(jdimn)=dist(jdimn)+(coord(jdimn,jpoin)-center(jdimn))/nnode
                    end do
                end do
            endif

            disgi=0.
            do idimn=1,ndimn
                disgi(idimn,idimn)=1.
            end do

            if(ndimn==2)then
                disgi(1,3)=-dist(2)
                disgi(2,3)= dist(1)
            elseif(ndimn==3)then
                disgi(1,5)= dist(3)
                disgi(1,6)=-dist(2)

                disgi(2,4)=-dist(3)
                disgi(2,6)= dist(1)

                disgi(3,4)= dist(2)
                disgi(3,5)=-dist(1)
            endif

            disli=0.
            if(block_stab/=1)then
                disli=rot.x.disgi
                disli=disli*coef
            elseif(block_stab==1)then
                disli(1:ndimn,:)=rot.x.disgi
                if(ndimn==2)disli(3,3)=1.
                if(ndimn==3) &
                    disli(ndimn+1:3*(ndimn-1),ndimn+1:3*(ndimn-1))=rot  !2015/11/17
                disli=disli*coef

                !write(7,*)'igapb=',igapb,'kpoin=',kpoin,'rot=',rot,'coef=',coef,'disli=',disli
            endif




            do idimn=1,gapb(igapb)%nrdof
                itotvbt=npgblock*kdimn+idimn
                do jdimn=1,kdimn
                    jtotvbt=(kpoin-1)*kdimn+jdimn
                    gapb(igapb)%cmatrix(itotvbt,jtotvbt)=disli(jdimn,idimn)   !2010/5/3
                    !write(7,*)'itotvbt=',itotvbt,'jtotvbt=',jtotvbt,'cmatrix=',gapb(igapb)%cmatrix(itotvbt,jtotvbt)
                enddo
            enddo !idimn

        enddo !kpoin




45      format(15e15.3)

        if(type_problem=='F')then
            fact=1+damp_ctt*theta1*ditime
            do idimn=1,gapb(igapb)%nrdof
                itotvbt=npgblock*ndimn+idimn
                gapb(igapb)%cmatrix(itotvbt,itotvbt)=-mass_inertia(idimn)*fact
                gapb(igapb)%rstiff(idimn,idimn)=-mass_inertia(idimn) !20121216
            enddo !idimn
        endif

        !nullify(listrdof)
    enddo !igapb
    deallocate(rot,disgi,disli,dist,center)
    if(type_problem=='F')deallocate(mass_inertia)

    end subroutine forAdirect_back_analysis   !20150925

    subroutine forAdirect !fzx !形成A矩阵  2010/7/13

    integer(ink) igapb,npgblock,onetwo,ipoin,jpoin,ipair,idimn,itotvbt,itotv,ielem, &
        kpoin,lpoin,jdimn,jtotv,jtotvbt,nevab,ieqx,ievab,igaps,nnode,ii,matno,index,order_int,ngaus, &
        ielgroup,igaus,igroup,jgroup
    real   (irk) coef,tvol,density,djacb,rr,fact,thick
    real   (irk),allocatable::eldis(:),eload(:),rot(:,:),disgi(:,:),disli(:,:),dist(:),center(:),mass_inertia(:),gpcod(:)
    integer(ink),pointer::ldofs(:)
    real   (irk),pointer::estif(:,:)


    allocate(center(ndimn),rot(ndimn,ndimn),dist(ndimn),disgi(ndimn,3*(ndimn-1)))
    if(block_stab/=1) allocate(disli(ndimn,3*(ndimn-1)))
    if(block_stab==1) allocate(disli(3*(ndimn-1),3*(ndimn-1)))
    if(type_problem=='F')allocate(mass_inertia(3*(ndimn-1)))
    !		   write(7,*)'npdisp in formatrixa *************jpoin,ipoin,gapb(igapb)%npdisp(:,jpoin,idimn)'

    write(7,*)'foradirect'

    allocate(gpcod(ndimn))
    do igapb=1,ngapb
        if(gapb(igapb)%nrdof==0) cycle   !tcl
        gapb(igapb)%center=0.;
        if(type_problem=='F')gapb(igapb)%mass_inertia=0.
        tvol=0.
        do igroup=1,gapb(igapb)%ngroupb
            jgroup=gapb(igapb)%listgroupb(igroup)
            if(jgroup<0)cycle  !20191031
            matno = group(jgroup)%matno
            index = group(jgroup)%index

            thick=1.
            if (ndimn==2.or.index==22.or.index==26)thick  =props(matno)%mechanical%solid%thickness     !20230910
            density=thick*props(matno)%mechanical%solid%density
            order_int=elkn(index)%el_field(1)%order_intrules(1)

            ngaus = elkn(index)%ggaus(order_int)%ngaus
            do ielgroup=1,group(jgroup)%nelgroup
                ielem = group(jgroup)%list(ielgroup)

                do igaus=1,ngaus
                    djacb=element(ielem)%egaus(order_int)%djacb(igaus)
                    gpcod=element(ielem)%egaus(order_int)%gpcod(:,igaus)
                    gapb(igapb)%center=gapb(igapb)%center+djacb*gpcod*density
                    tvol=tvol+djacb*density
                end do
            enddo
        end do
        gapb(igapb)%center=gapb(igapb)%center/tvol
        write(7,*)'tvol=',tvol,'igapb=',igapb,'center=', gapb(igapb)%center

        if(type_problem=='F')then
            do igroup=1,gapb(igapb)%ngroupb
                jgroup=gapb(igapb)%listgroupb(igroup)
                if(jgroup<0)cycle  !20191031
                matno = group(jgroup)%matno
                index = group(jgroup)%index
                thick=1.
                if (ndimn==2.or.index==22.or.index==26)thick  =props(matno)%mechanical%solid%thickness
                density=thick*props(matno)%mechanical%solid%density
                order_int=elkn(index)%el_field(1)%order_intrules(1)
                ngaus = elkn(index)%ggaus(order_int)%ngaus
                do ielgroup=1,group(jgroup)%nelgroup
                    ielem = group(jgroup)%list(ielgroup)

                    do igaus=1,ngaus
                        djacb=element(ielem)%egaus(order_int)%djacb(igaus)
                        gpcod=element(ielem)%egaus(order_int)%gpcod(:,igaus)
                        if(ndimn==2)then
                            rr=(gpcod(1)-gapb(igapb)%center(1))**2+(gpcod(2)-gapb(igapb)%center(2))**2
                            gapb(igapb)%mass_inertia(1)=gapb(igapb)%mass_inertia(1)+djacb*density
                            gapb(igapb)%mass_inertia(2)=gapb(igapb)%mass_inertia(2)+djacb*density
                            gapb(igapb)%mass_inertia(3)=gapb(igapb)%mass_inertia(3)+djacb*rr*density
                        elseif(ndimn==3)then
                            gapb(igapb)%mass_inertia(1)=gapb(igapb)%mass_inertia(1)+djacb*density
                            gapb(igapb)%mass_inertia(2)=gapb(igapb)%mass_inertia(2)+djacb*density
                            gapb(igapb)%mass_inertia(3)=gapb(igapb)%mass_inertia(3)+djacb*density
                            rr=(gpcod(2)-gapb(igapb)%center(2))**2+(gpcod(3)-gapb(igapb)%center(3))**2
                            gapb(igapb)%mass_inertia(4)=gapb(igapb)%mass_inertia(4)+djacb*rr*density
                            rr=(gpcod(1)-gapb(igapb)%center(1))**2+(gpcod(3)-gapb(igapb)%center(3))**2
                            gapb(igapb)%mass_inertia(5)=gapb(igapb)%mass_inertia(5)+djacb*rr*density
                            rr=(gpcod(2)-gapb(igapb)%center(2))**2+(gpcod(1)-gapb(igapb)%center(1))**2
                            gapb(igapb)%mass_inertia(6)=gapb(igapb)%mass_inertia(6)+djacb*rr*density
                        endif

                    end do
                enddo
            end do
            write(7,*)'igapb=',igapb,'gapb(igapb)%mass_inertia=', gapb(igapb)%mass_inertia
        endif
    enddo

    deallocate(gpcod)


    do igapb=1,ngapb

        if(gapb(igapb)%nrdof==0) cycle   !tcl
        npgblock=gapb(igapb)%npgblock
        center=gapb(igapb)%center
        if(type_problem=='F')mass_inertia=gapb(igapb)%mass_inertia

        !write(7,*)'npdisp++++++++++++++'
        do jpoin=1,gapb(igapb)%npblock
            ipoin=gapb(igapb)%nodeblock(jpoin)
            dist=0.
            do jdimn=1,ndimn
                dist(jdimn)=coord(jdimn,ipoin)-center(jdimn)
            end do
            disgi=0.
            do idimn=1,ndimn
                disgi(idimn,idimn)=1.
            end do

            if(ndimn==2)then
                disgi(1,3)=-dist(2)
                disgi(2,3)= dist(1)
            elseif(ndimn==3)then
                disgi(1,5)= dist(3)
                disgi(1,6)=-dist(2)

                disgi(2,4)=-dist(3)
                disgi(2,6)= dist(1)

                disgi(3,4)= dist(2)
                disgi(3,5)=-dist(1)
            endif

            do idimn=1,gapb(igapb)%nrdof
                do jdimn=1,ndimn
                    gapb(igapb)%npdisp(jdimn,jpoin,idimn)=disgi(jdimn,idimn) !存每个点的位移，为计算A(T)F
                enddo
            end do
            if(block_stab==1)then
                do jdimn=ndimn+1,3*(ndimn-1)
                    gapb(igapb)%npdisp(jdimn,jpoin,jdimn)=1. !存每个点的位移，为计算A(T)F
                enddo
            endif
            !		  write(7,*)jpoin,ipoin,gapb(igapb)%npdisp(:,jpoin,:)
        enddo


        !do idimn=1,gapb(igapb)%nrdof
        !         gapb(igapb)%rstiff(idimn,idimn)=1.e-5 !20121216
        !enddo !idimn


        ! if(type_problem=='F')then
        ! do idimn=1,gapb(igapb)%nrdof
        !             gapb(igapb)%rstiff(idimn,idimn)=-mass_inertia(idimn) !20121216
        !    enddo !idimn
        !endif

        if (restart_ctt/=0) cycle

        do kpoin=1,npgblock
            onetwo=gapb(igapb)%nodegblock_onetwo(kpoin) !前四个点还是后四个点，也就是第一点还是第二点
            if(onetwo==1)coef=1.
            if(onetwo==2)coef=-1.
            igaps=gapb(igapb)%nodegblock_igaps(kpoin)
            ipair=gapb(igapb)%nodegblock_ipairs(kpoin)
            rot=gaps(igaps)%rot(:,:,ipair)


            dist=0.
            if(contactpe==1)then
                lpoin=gapb(igapb)%nodegblock(kpoin)
                do jdimn=1,ndimn
                    dist(jdimn)=coord(jdimn,lpoin)-center(jdimn)
                end do
            elseif(contactpe==2)then
                nnode=2
                if(ndimn==3)nnode=4
                do ii=1,nnode
                    jpoin=gaps(igaps)%pairnode((onetwo-1)*nnode+ii,ipair)
                    do jdimn=1,ndimn
                        dist(jdimn)=dist(jdimn)+(coord(jdimn,jpoin)-center(jdimn))/nnode
                    end do
                end do
            endif

            disgi=0.
            do idimn=1,ndimn
                disgi(idimn,idimn)=1.
            end do

            if(ndimn==2)then
                disgi(1,3)=-dist(2)
                disgi(2,3)= dist(1)
            elseif(ndimn==3)then
                disgi(1,5)= dist(3)
                disgi(1,6)=-dist(2)

                disgi(2,4)=-dist(3)
                disgi(2,6)= dist(1)

                disgi(3,4)= dist(2)
                disgi(3,5)=-dist(1)
            endif

            disli=0.
            if(block_stab/=1)then
                disli=rot.x.disgi
                disli=disli*coef
            elseif(block_stab==1)then
                disli(1:ndimn,:)=rot.x.disgi
                if(ndimn==2)disli(3,3)=1.
                if(ndimn==3) &
                    disli(ndimn+1:3*(ndimn-1),ndimn+1:3*(ndimn-1))=rot
                disli=disli*coef

                !write(7,*)'igapb=',igapb,'kpoin=',kpoin,'rot=',rot,'coef=',coef,'disli=',disli
            endif



            if(block_stab/=1)then
                do idimn=1,gapb(igapb)%nrdof
                    itotvbt=npgblock*ndimn+idimn
                    do jdimn=1,ndimn
                        jtotvbt=(kpoin-1)*ndimn+jdimn
                        gapb(igapb)%cmatrix(jtotvbt,itotvbt)=disli(jdimn,idimn)
                        gapb(igapb)%cmatrix(itotvbt,jtotvbt)=disli(jdimn,idimn)   !2010/5/3
                    enddo
                enddo !idimn
            elseif(block_stab==1)then
                do idimn=1,gapb(igapb)%nrdof
                    itotvbt=npgblock*3*(ndimn-1)+idimn
                    do jdimn=1,3*(ndimn-1)
                        jtotvbt=(kpoin-1)*3*(ndimn-1) +jdimn
                        gapb(igapb)%cmatrix(jtotvbt,itotvbt)=disli(jdimn,idimn)
                        gapb(igapb)%cmatrix(itotvbt,jtotvbt)=disli(jdimn,idimn)   !2010/5/3
                    enddo
                enddo !idimn
            endif

        enddo !kpoin

        !write(7,*)'igapb=',igapb
        !   write(7,*)'cmatrix='
        !   do itotvbt=1,gapb(igapb)%ntotv_bt
        !       write(7,45)gapb(igapb)%cmatrix(itotvbt,:)
        !   end do
45      format(15e15.3)

        if(type_problem=='F'.and.gapb(igapb)%eblock==0)then
            fact=1+damp_ctt*theta1*ditime
            do idimn=1,gapb(igapb)%nrdof
                itotvbt=npgblock*ndimn+idimn
                gapb(igapb)%cmatrix(itotvbt,itotvbt)=-mass_inertia(idimn)*fact
                gapb(igapb)%rstiff(idimn,idimn)=-mass_inertia(idimn) !20121216
            enddo !idimn
        endif


        !nullify(listrdof)
    enddo !igapb
    deallocate(rot,disgi,disli,dist,center)
    if(type_problem=='F')deallocate(mass_inertia)

    end subroutine forAdirect   !tcl

    subroutine matrix_rigid_dis   !20211121
    character(20) text
    integer(ink) igapb,npgblock,ipoin,jpoin,idimn,itotvbt,itotv,ielem, &
        kpoin,lpoin,jdimn,jtotv,jtotvbt,nevab,ieqx,ievab,igaps,nnode,ii,matno,index,order_int,ngaus, &
        ielgroup,igaus,igroup,jgroup,igdis_bk,ngroup_bk,npoin_bk
    integer(ink),allocatable:: appear_gdis_bk(:)
    integer(ink),pointer:: lnods(:)
    real   (irk) tvol,density,djacb,fact,thick
    real   (irk),allocatable::disgi(:,:),dist(:),center(:),gpcod(:),cor(:)

    read(back_ctl_unit,*)text
    print *,'text=',text
    read(back_ctl_unit,*)ngdis_bk
    print *,'ngdis_bk=',ngdis_bk
    allocate(rigid_bk(ngdis_bk),appear_gdis_bk(ngroup))

    read(back_ctl_unit,*)text
    do igdis_bk=1,ngdis_bk
        allocate(rigid_bk(igdis_bk)%fixed_dis(3*(ndimn-1)))
        read(back_ctl_unit,*)rigid_bk(igdis_bk)%fixed_dis
        print *,'fixed_dis=',rigid_bk(igdis_bk)%fixed_dis
        read(back_ctl_unit,*)appear_gdis_bk(:)
        print *,'appear_gdis_bk=',appear_gdis_bk
        read(back_ctl_unit,*)npoin_bk
        rigid_bk(igdis_bk)%npoin_bk=npoin_bk
        allocate(rigid_bk(igdis_bk)%node_bk(npoin_bk))
        read(back_ctl_unit,*)rigid_bk(igdis_bk)%node_bk  !对应的是para_points(1:npoints_pb)的序号
        ngroup_bk=sum(appear_gdis_bk)
        print *,'ngroup_bk=',ngroup_bk
        rigid_bk(igdis_bk)%ngroup_bk=ngroup_bk
        allocate(rigid_bk(igdis_bk)%group_bk(ngroup_bk))
        ngroup_bk=0
        do igroup=1,ngroup
            if(appear_gdis_bk(igroup)==0)cycle
            ngroup_bk=ngroup_bk+1
            rigid_bk(igdis_bk)%group_bk(ngroup_bk)=igroup
        end do
        print *,'ngroup_bk=',ngroup_bk
        print *,'rigid_bk(igdis_bk)%group_bk=',rigid_bk(igdis_bk)%group_bk

    end do

    allocate(gpcod(ndimn),cor(ndimn),center(ndimn),dist(ndimn),disgi(ndimn,3*(ndimn-1)))

    do igdis_bk=1,ngdis_bk

        allocate(rigid_bk(igdis_bk)%center(ndimn))
        rigid_bk(igdis_bk)%center=0.
        tvol=0.
        do jgroup=1,rigid_bk(igdis_bk)%ngroup_bk
            igroup=rigid_bk(igdis_bk)%group_bk(jgroup)
            !print *,'jgroup=',igroup
            matno = group(igroup)%matno
            index = group(igroup)%index

            thick=1.
            if (ndimn==2.or.index==22.or.index==26)thick  =props(matno)%mechanical%solid%thickness
            density=thick*props(matno)%mechanical%solid%density
            order_int=elkn(index)%el_field(1)%order_intrules(1)
            ngaus = elkn(index)%ggaus(order_int)%ngaus
            do ielgroup=1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                lnods=>element(ielem)%field(1)%lnods_f

                nullify(lnods)
                do igaus=1,ngaus
                    djacb=element(ielem)%egaus(order_int)%djacb(igaus)
                    gpcod=element(ielem)%egaus(order_int)%gpcod(:,igaus)
                    rigid_bk(igdis_bk)%center=rigid_bk(igdis_bk)%center+djacb*gpcod*density
                    tvol=tvol+djacb*density
                end do
            enddo  !ielgroup
        end do  !jgroup
        rigid_bk(igdis_bk)%center=rigid_bk(igdis_bk)%center/tvol
    enddo  !igdis_bk
    deallocate(gpcod)


    do igdis_bk=1,ngdis_bk
        npoin_bk=rigid_bk(igdis_bk)%npoin_bk
        allocate(rigid_bk(igdis_bk)%npdisp(ndimn,npoin_bk,3*(ndimn-1)))
        rigid_bk(igdis_bk)%npdisp=0.

        center=rigid_bk(igdis_bk)%center
        write(7,*)'igdis_bk=',igdis_bk,'center=',center

        !write(7,*)'npdisp++++++++++++++'
        do jpoin=1,npoin_bk
            ipoin=rigid_bk(igdis_bk)%node_bk(jpoin)
            cor=para_points(jpoin)%cor
            !write(7,*)'jpoin=',jpoin,'cor=',cor
            dist=0.
            do jdimn=1,ndimn
                dist(jdimn)=cor(jdimn)-center(jdimn)
            end do
            disgi=0.
            do idimn=1,ndimn
                disgi(idimn,idimn)=1.
            end do

            if(ndimn==2)then
                disgi(1,3)=-dist(2)
                disgi(2,3)= dist(1)
            elseif(ndimn==3)then
                disgi(1,5)= dist(3)
                disgi(1,6)=-dist(2)

                disgi(2,4)=-dist(3)
                disgi(2,6)= dist(1)

                disgi(3,4)= dist(2)
                disgi(3,5)=-dist(1)
            endif

            do idimn=1,3*(ndimn-1)
                do jdimn=1,ndimn
                    rigid_bk(igdis_bk)%npdisp(jdimn,jpoin,idimn)=disgi(jdimn,idimn) !存每个点的位移，为计算A(T)F
                enddo
            end do
            !		  write(7,*)jpoin,ipoin,gapb(igapb)%npdisp(:,jpoin,:)
        enddo

    enddo !igapb
    deallocate(disgi,dist,cor,center,appear_gdis_bk)
45  format(15e15.3)


    end subroutine matrix_rigid_dis   !20211121

    subroutine matrix_nodal_value   !20211201
    character(20) text
    integer(ink) ipoin,igroup,jgroup,ielem,ielgroup,igdis_bk,ngroup_bk, &
        npoin_bk,ndofn_bk
    integer(ink),allocatable:: appear_nodvar_bk(:),listnode(:)
    integer(ink),pointer:: lnods(:)

    read(back_ctl_unit,*)text
    read(back_ctl_unit,*)ngval_bk
    allocate(nodvar_bk(ngval_bk),appear_nodvar_bk(ngroup))

    read(back_ctl_unit,*)text
    do igdis_bk=1,ngval_bk
        read(back_ctl_unit,*)ndofn_bk
        nodvar_bk(igdis_bk)%ndofn_bk=ndofn_bk
        allocate(nodvar_bk(igdis_bk)%dof_bk(ndofn_bk))
        read(back_ctl_unit,*)nodvar_bk(igdis_bk)%dof_bk

        read(back_ctl_unit,*)appear_nodvar_bk(:)
        ngroup_bk=sum(appear_nodvar_bk)
        nodvar_bk(igdis_bk)%ngroup_bk=ngroup_bk
        allocate(nodvar_bk(igdis_bk)%group_bk(ngroup_bk))
        ngroup_bk=0
        do igroup=1,ngroup
            if(appear_nodvar_bk(igroup)==0)cycle
            ngroup_bk=ngroup_bk+1
            nodvar_bk(igdis_bk)%group_bk(ngroup_bk)=igroup
        end do
        print *,'ngroup_bk=',ngroup_bk
        print *,'nodvar_bk(igdis_bk)%group_bk=',nodvar_bk(igdis_bk)%group_bk

    end do

    allocate (listnode(npoin))
    do igdis_bk=1,ngval_bk
        listnode=0

        do jgroup=1,nodvar_bk(igdis_bk)%ngroup_bk
            igroup=nodvar_bk(igdis_bk)%group_bk(jgroup)
            do ielgroup=1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                lnods=>element(ielem)%field(1)%lnods_f
                listnode(lnods)=1
                nullify(lnods)
            enddo  !ielgroup
        end do  !jgroup
        npoin_bk=sum(listnode)
        nodvar_bk(igdis_bk)%npoin_bk=npoin_bk
        allocate(nodvar_bk(igdis_bk)%node_bk(npoin_bk),nodvar_bk(igdis_bk)%nodet_bk(npoin))

        npoin_bk=0
        nodvar_bk(igdis_bk)%nodet_bk=0
        do ipoin=1,npoin
            if(listnode(ipoin)==0)cycle
            npoin_bk=npoin_bk+1
            nodvar_bk(igdis_bk)%node_bk(npoin_bk)=ipoin
            nodvar_bk(igdis_bk)%nodet_bk(ipoin)=npoin_bk
        end do
    enddo  !igdis_bk

    deallocate(appear_nodvar_bk,listnode)


    end subroutine matrix_nodal_value   !20211201
