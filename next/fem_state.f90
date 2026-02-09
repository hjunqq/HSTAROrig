    SUBROUTINE gpvarupdate

    character(10) fieldid,class,name,material,model
    integer(ink) igroup,ielgroup,ielem,matno,order_int,index,icr,icreep  !20180630
    DO igroup =1,ngroup   ! --1
        if (appear(igroup)>0) then
            ! get information from the group level
            fieldid=group(igroup)%fieldid
            class  =group(igroup)%class
            !         if (fieldid(1:1)=='U'.and.class=='CO')then
            if (fieldid(1:1)=='U')then !--2
                matno = group(igroup)%matno
                index = group(igroup)%index
                name  =props(matno)%name
                material=props(matno)%mechanical%solid%material
                if(material=='GOODMAN') &
                    model=props(matno)%mechanical%solid%Goodman%model


                icreep =props(matno)%mechanical%solid%icreep  !20180630
                icr=0
                if(material=='CONCRETE')icr=props(matno)%mechanical%solid%Concrete%icr
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    if(ice0(ielem)==1) goto 100
                    element(ielem)%field(1)%gpvar0=element(ielem)%field(1)%gpvar
                    if(index/=20.and.index/=21)  &   !20211125
                        element(ielem)%field(1)%sigz=element(ielem)%field(1)%gpvar(ndimn,:) !ep2010
                    if(icr==2.or.icr==3.or.icr==5.or.icr==6)  & !zhao09
                        element(ielem)%field(1)%strain0=element(ielem)%field(1)%strain

                    if(material=='GOODMAN'.and.model=='FCM')  & !20210126
                        element(ielem)%field(1)%strain0=element(ielem)%field(1)%strain

                    if(icreep==4)  & !20180630
                        element(ielem)%field(1)%vkstrain0=element(ielem)%field(1)%vkstrain

                    !! contact
                    if (name=='CONTACT') then
                        element(ielem)%field(1)%gapg0=element(ielem)%field(1)%gapg
                        element(ielem)%field(1)%gapn0=element(ielem)%field(1)%gapn
                        element(ielem)%field(1)%state0=element(ielem)%field(1)%state
                    elseif(name=='CRACK')then !crack 2006
                        element(ielem)%field(1)%state0=element(ielem)%field(1)%state
                    endif
                    !20231215_YL
                    !!20231007 止水
                    if(material=='GOODMAN') then
                        model=props(matno)%mechanical%solid%Goodman%model
                        if(model=='WATERTIGHT')then
                            element(ielem)%field(1)%relat_dis_nod0 =element(ielem)%field(1)%relat_dis_nod
                            element(ielem)%field(1)%relat_dis_gaus0=element(ielem)%field(1)%relat_dis_gaus
                        endif
                    endif
                    !20231215_YL


                    !!contact
100                 continue
                end do
            endif




            if (fieldid(1:2)=='UW') then !--3

                if(material=='SandPZ'.or.material=='ClayPZ'.or.material=='SoilPZ')then
                    order_int=elkn(index)%el_field(1)%order_intrules(1)
                    DO ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        if (ice0(ielem)/=1)then
                            element(ielem)%egaus(order_int)%iload0=element(ielem)%egaus(order_int)%iload
                            element(ielem)%egaus(order_int)%vdval0=element(ielem)%egaus(order_int)%vdval
                        endif
                    end do
                end if
            endif !--3
        endif !--2
    end do !--1

    END  SUBROUTINE gpvarupdate

    SUBROUTINE gpvarupdate1
    character(10) fieldid,class,name,material
    integer(ink) igroup,ielgroup,ielem,matno,order_int,index,icr
    DO igroup =1,ngroup
        if (appear(igroup)>0) then
            ! get information from the group level
            fieldid=group(igroup)%fieldid
            class  =group(igroup)%class
            !      if(fieldid(1:1)=='U'.and.class=='CO')then
            if (fieldid(1:1)=='U')then
                matno = group(igroup)%matno
                index = group(igroup)%index
                name    =props(matno)%name
                material=props(matno)%mechanical%solid%material
                icr=0
                if (material=='CONCRETE')icr=props(matno)%mechanical%solid%Concrete%icr
                DO ielgroup = 1,group1(igroup)%nelgroup
                    ielem = group1(igroup)%list(ielgroup)
                    if (jce1(ielem)/=1)then
                        element1(ielem)%field(1)%gpvar0=element1(ielem)%field(1)%gpvar
                        if (icr==2.or.icr==3.or.icr==5)  & !zhao09
                            element1(ielem)%field(1)%strain0=element1(ielem)%field(1)%strain

                        !! contact
                        if (name=='CONTACT') then
                            element1(ielem)%field(1)%gapg0=element1(ielem)%field(1)%gapg
                            element1(ielem)%field(1)%gapn0=element1(ielem)%field(1)%gapn
                            element1(ielem)%field(1)%state0=element1(ielem)%field(1)%state
                        endif
                        !!contact
                    endif
                end do
            endif


            if (fieldid(1:2)=='UW') then
                if(material=='SandPZ'.or.material=='ClayPZ'.or.material=='SoilPZ')then
                    order_int=elkn(index)%el_field(1)%order_intrules(1)
                    DO ielgroup = 1,group1(igroup)%nelgroup
                        ielem = group1(igroup)%list(ielgroup)
                        if (jce1(ielem)==1)then
                            element1(ielem)%egaus(order_int)%iload0=      &
                                element1(ielem)%egaus(order_int)%iload
                            element1(ielem)%egaus(order_int)%vdval0=      &
                                element1(ielem)%egaus(order_int)%vdval
                        endif
                    end do
                endif
            endif

        endif
    end do
    END  SUBROUTINE gpvarupdate1

    SUBROUTINE gpvarupdate2
    character(10) fieldid,class,name,material
    integer(ink) igroup,ielgroup,ielem,matno,order_int,index,icr
    DO igroup =1,ngroup
        if (appear(igroup)>0) then
            ! get information from the group level
            fieldid=group(igroup)%fieldid
            class  =group(igroup)%class
            !      if(fieldid(1:1)=='U'.and.class=='CO')then
            if (fieldid(1:1)=='U')then
                matno = group(igroup)%matno
                index = group(igroup)%index
                name    =props(matno)%name
                material=props(matno)%mechanical%solid%material
                icr=0
                if (material=='CONCRETE')icr=props(matno)%mechanical%solid%Concrete%icr
                DO ielgroup = 1,group2(igroup)%nelgroup
                    ielem = group2(igroup)%list(ielgroup)
                    element2(ielem)%field(1)%gpvar0=element2(ielem)%field(1)%gpvar
                    if (icr==2.or.icr==3.or.icr==5)  & !zhao09
                        element2(ielem)%field(1)%strain0=element2(ielem)%field(1)%strain

                    !! contact
                    if (name=='CONTACT') then
                        element2(ielem)%field(1)%gapg0=element2(ielem)%field(1)%gapg
                        element2(ielem)%field(1)%gapn0=element2(ielem)%field(1)%gapn
                        element2(ielem)%field(1)%state0=element2(ielem)%field(1)%state
                    endif
                    !!contact
                end do
            endif


            if (fieldid(1:2)=='UW') then
                if(material=='SandPZ'.or.material=='ClayPZ'.or.material=='SoilPZ')then
                    order_int=elkn(index)%el_field(1)%order_intrules(1)
                    DO ielgroup = 1,group2(igroup)%nelgroup
                        ielem = group2(igroup)%list(ielgroup)
                        element2(ielem)%egaus(order_int)%iload0=      &
                            element2(ielem)%egaus(order_int)%iload
                        element2(ielem)%egaus(order_int)%vdval0=      &
                            element2(ielem)%egaus(order_int)%vdval
                    end do
                endif
            endif

        endif
    end do
    END  SUBROUTINE gpvarupdate2

    SUBROUTINE gpvar_initial
    character(10) fieldid,class,name,material
    integer(ink) igroup,ielgroup,ielem,matno,order_int,index
    DO igroup =1,ngroup
        !if (appear(igroup)==-1) then
        ! get information from the group level
        fieldid=group(igroup)%fieldid
        class  =group(igroup)%class
        !      if(fieldid(1:1)=='U'.and.class=='CO')then
        if (fieldid(1:1)=='U')then
            DO ielgroup = 1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                element(ielem)%field(1)%gpvar0=0.0
                element(ielem)%field(1)%gpvar =0.0
                element(ielem)%stres0=0.
                if (associated(element(ielem)%field(1)%strain0)) &
                    element(ielem)%field(1)%strain0=0.
                if (associated(element(ielem)%field(1)%strain )) &
                    element(ielem)%field(1)%strain =0.
            end do
        endif
        if (fieldid(1:2)=='UW') then
            matno = group(igroup)%matno
            index = group(igroup)%index
            material=props(matno)%mechanical%solid%material
            if(material=='SandPZ'.or.material=='ClayPZ'.or.material=='SoilPZ')then
                order_int=elkn(index)%el_field(1)%order_intrules(1)
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    element(ielem)%egaus(order_int)%iload0=0.
                    element(ielem)%egaus(order_int)%vdval0=0.
                    element(ielem)%egaus(order_int)%iload =0.
                    element(ielem)%egaus(order_int)%vdval =0.
                end do
            endif
        endif

        !endif
    end do
    END  SUBROUTINE gpvar_initial

    SUBROUTINE gpvar1_initial
    character(10) fieldid,class,name,material
    integer(ink) igroup,ielgroup,ielem,matno,order_int,index
    DO igroup =1,ngroup
        if (appear(igroup)==-1) then
            ! get information from the group level
            fieldid=group(igroup)%fieldid
            class  =group(igroup)%class
            !      if(fieldid(1:1)=='U'.and.class=='CO')then
            if (fieldid(1:1)=='U')then
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    element(ielem)%field(1)%gpvar0=0.0
                    element(ielem)%field(1)%gpvar =0.0
                    element(ielem)%stres0=0.
                    if (associated(element(ielem)%field(1)%strain0)) &
                        element(ielem)%field(1)%strain0=0.
                    if (associated(element(ielem)%field(1)%strain )) &
                        element(ielem)%field(1)%strain =0.
                end do
            endif
            if (fieldid(1:2)=='UW') then
                matno = group(igroup)%matno
                index = group(igroup)%index
                material=props(matno)%mechanical%solid%material
                if(material=='SandPZ'.or.material=='ClayPZ'.or.material=='SoilPZ')then
                    order_int=elkn(index)%el_field(1)%order_intrules(1)
                    DO ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        element(ielem)%egaus(order_int)%iload0=0.
                        element(ielem)%egaus(order_int)%vdval0=0.
                        element(ielem)%egaus(order_int)%iload =0.
                        element(ielem)%egaus(order_int)%vdval =0.
                    end do
                endif
            endif

        endif
    end do
    END  SUBROUTINE gpvar1_initial

    SUBROUTINE gpvar2_initial
    character(10) fieldid,class,name
    integer(ink) igroup,ielgroup,ielem,matno,order_int,index
    DO igroup =1,ngroup
        if (appear(igroup)==2) then
            ! get information from the group level
            fieldid=group(igroup)%fieldid
            class  =group(igroup)%class
            !      if(fieldid(1:1)=='U'.and.class=='CO')then
            if (fieldid(1:1)=='U')then
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    element(ielem)%field(1)%gpvar0=0.0
                    element(ielem)%field(1)%gpvar =0.0
                    element(ielem)%stres0=0.
                    if (associated(element(ielem)%field(1)%strain0)) &
                        element(ielem)%field(1)%strain0=0.
                    if (associated(element(ielem)%field(1)%strain )) &
                        element(ielem)%field(1)%strain =0.
                end do
            endif
            if (fieldid(1:2)=='UW') then
                matno = group(igroup)%matno
                index = group(igroup)%index
                material=props(matno)%mechanical%solid%material
                if(material=='SandPZ'.or.material=='ClayPZ'.or.material=='SoilPZ')then
                    order_int=elkn(index)%el_field(1)%order_intrules(1)
                    DO ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        element(ielem)%egaus(order_int)%iload0=0.
                        element(ielem)%egaus(order_int)%vdval0=0.
                        element(ielem)%egaus(order_int)%iload =0.
                        element(ielem)%egaus(order_int)%vdval =0.
                    end do
                endif
            endif

        endif
    end do
    END  SUBROUTINE gpvar2_initial

    subroutine conver_load

    real   (irk),allocatable:: refor(:),toforx(:),refory(:),tofory(:)
    integer(ink) igroup,nrfields,nelgroup,ielgroup,ielem,imcon,                &
        ifield,first_node, second_node,npairs,ipair,            &
        itotv,jtotv,order_freedom,ilink,ndofn,                  &
        inode,ipoin,np_unode,ii,iedge,idimn,felem,i0, &
        i1,ipairs,npairs_wc,nline_g_w,twater_curve,iwc,  &
        igapbf,igapb,npgblock,ij,igaps,kkdimn

    integer(ink) iintf,nintf
    character(10)fieldid
    integer(ink),pointer::ldofs(:),link_freedom(:),pairnode(:,:),lnods(:),pairnode_wc(:)
    real   (irk),pointer::eload(:),value(:),ks(:,:)
    real   (irk) tfi,resid,retot,max_refor,max_tofor,ratio1,ratio2,xxxx,coef1,dfact
    real   (irk),allocatable::heat_wc(:)  !20210417


    !   Checks convergence for load part
    print *,'                    istep=',istep,'iiter=',iiter
    print *,'                    conver_load check'
    write(chkunit,*)
    write(chkunit,*)'        Convergence check for Residual force,ttime=,',ttime
    allocate(refor(ntotv),toforx(ntotv),refory(ntotv),tofory(ntotv))
    nchek=0
    resid=0.0
    retot=0.0
    stfor=0.0
    toforx=0.0
    refor=0.0
    tofory=0.0
    refory=0.0


    do igroup=1,ngroup
        if (appear(igroup)>0) then   !20191006
            ! get information from the group level
            nrfields =group(igroup)%nrfields
            nelgroup =group(igroup)%nelgroup
            do ielgroup=1,nelgroup
                ielem=group(igroup)%list(ielgroup)
                if (ice0(ielem)/=1)then
                    do ifield=1,nrfields
                        if (associated(element(ielem)%field(ifield)%eload)) then
                            ldofs=> element(ielem)%field(ifield)%ldofs_f
                            eload=> element(ielem)%field(ifield)%eload

                            ndofn=size(ldofs)
                            do itotv=1,ndofn

                                stfor(ldofs(itotv))=stfor(ldofs(itotv))+eload(itotv)
                                ! if(ldofs(itotv)==5218) &
                                !write(7,*)'ielem=',ielem,'ifield=',ifield,'idofn=',itotv,'stfor=',stfor(ldofs(itotv)),'eload=',eload(itotv)
                            end do
                            nullify(ldofs,eload)
                        endif
                    end do   !! for ifield
                endif
            end do  !! for ielgroup
        end if   !! for appear
    end do !! for igroup
    !!!!!
    write(7,*)'total force on spring nodes of dam_foundation interface'

    !ifs2006 zhao, 06/03/29 , icaddmass

    if (icaddmass/=0)then
        do ipoin=1,npoin
            if(icmp(ipoin)==0)cycle
            do idimn=1,ndimn
                itotv=nodfn(idimn,ipoin)
                if(itotv==0)cycle
                xxxx=result_second(itotv)
                stfor(itotv)=stfor(itotv)+addmp(idimn,ipoin)*xxxx
            enddo
        enddo
    endif
    !2013/4/12

    if (nmcon/=0)then
        do imcon=1,nmcon
            ipoin=lmcon(imcon)
            do idimn=1,ndimn
                itotv=nodfn(idimn,ipoin)
                if(itotv==0)cycle
                xxxx=result_second(itotv)
                stfor(itotv)=stfor(itotv)+rmcon(idimn,imcon)*xxxx
            enddo
        enddo
    endif

    if (nbspring>0)then  !20150925
        do imcon=1,nbspring
            itotv=bspring(imcon)%listdof
            if(itotv==0)cycle
            stfor(itotv)=stfor(itotv)+bspring(imcon)%eload
        enddo
    endif      !20150925

    !!!!1
    if (rmesh>0.and.nelem1>0)then
        do igroup=1,ngroup
            if (appear(igroup)>0) then
                ! get information from the group level
                nrfields =group(igroup)%nrfields
                nelgroup =group1(igroup)%nelgroup
                do ielgroup=1,nelgroup
                    ielem=group1(igroup)%list(ielgroup)
                    if (jce1(ielem)/=1)then
                        do ifield=1,nrfields
                            if (associated(element1(ielem)%field(ifield)%eload)) then
                                ldofs=> element1(ielem)%field(ifield)%ldofs_f
                                eload=> element1(ielem)%field(ifield)%eload

                                ndofn=size(ldofs)
                                do itotv=1,ndofn
                                    stfor(ldofs(itotv))=stfor(ldofs(itotv))+eload(itotv)

                                end do
                                nullify(ldofs,eload)
                            endif
                        end do   !! for ifield
                    endif
                end do  !! for ielgroup
            end if   !! for appear
        end do !! for igroup
    endif

    !!!!2
    if (rmesh>1.and.nelem2>0)then
        do igroup=1,ngroup
            if (appear(igroup)>0) then
                ! get information from the group level
                nrfields =group(igroup)%nrfields
                nelgroup =group2(igroup)%nelgroup
                do ielgroup=1,nelgroup
                    ielem=group2(igroup)%list(ielgroup)
                    do ifield=1,nrfields
                        if (associated(element2(ielem)%field(ifield)%eload)) then
                            ldofs=> element2(ielem)%field(ifield)%ldofs_f
                            eload=> element2(ielem)%field(ifield)%eload

                            ndofn=size(ldofs)
                            do itotv=1,ndofn
                                stfor(ldofs(itotv))=stfor(ldofs(itotv))+eload(itotv)
                            end do
                            nullify(ldofs,eload)
                        endif
                    end do   !! for ifield
                end do  !! for ielgroup
            end if   !! for appear
        end do !! for igroup
    endif

    !!ifs2000
    do ielem=1,nifsgroup
        ldofs=> tifs(ielem)%ldofs
        eload=> tifs(ielem)%eload
        stfor(ldofs)=stfor(ldofs)+eload
        nullify(ldofs,eload)
    end do

    do ielem=1,nabssgroup
        ldofs=> tabss(ielem)%ldofs
        eload=> tabss(ielem)%eload
        stfor(ldofs)=stfor(ldofs)+eload
        nullify(ldofs,eload)
    end do
    !!ifs2000

    !ifs2006 zhao, 06/03/29
    if(icaddmass==0)then    !20231215YL 对附加质量法不需要以下集成
        do iedge=1,ifsnedge
            felem=ifsedges(iedge)%felem
            igroup=element(felem)%group
            if(appear(igroup)==0)cycle
            ldofs=>ifsedges(iedge)%ldofs
            eload=>ifsedges(iedge)%eload
            stfor(ldofs)=stfor(ldofs)+eload
            nullify(ldofs,eload)
        enddo
    endif

    !! stablize
    if (stabpw==1) then
        do igroup=1,ngroup
            fieldid=group(igroup)%fieldid
            if (appear(igroup)>0.and.(fieldid(1:2)=='UP'.or.fieldid(1:2)=='UW')) then
                do ipoin=1,group(igroup)%np_unode
                    np_unode=group(igroup)%unode(ipoin)%np_unode
                    if (np_unode/=0) then
                        eload=>group(igroup)%unode(ipoin)%patch_load
                        lnods=>group(igroup)%unode(ipoin)%patch_nod
                        do inode=1,np_unode
                            itotv=nodfn(ndimn+1,lnods(inode))
                            stfor(itotv)=stfor(itotv)+eload(inode)
                        end do
                        nullify(eload,lnods)
                    endif
                end do
            endif
        end do
    endif
    !! end of stablize

    if (ground_inf/=0) then
        !  write(7,*)'itotv,ldofs_space,eload_space,stfor='
        !do itotv=1,size(ldofs_space)
        !    write(7,*)itotv,ldofs_space(itotv),eload_space(itotv),stfor(ldofs_space(itotv))
        !end do
        stfor(ldofs_space)=stfor(ldofs_space)+eload_space

    endif

    if(nwcpipe/=0)then !20210417
        !write(7,*)'nwcpipe:itotv,stfor,tofor='
        do i0=1,nwcpipe
            nline_g_w=wc_pipe(i0)%nline_g_w
            coef1= wc_pipe(i0)%iwc
            iwc=wc_pipe(i0)%iwc
            twater_curve=wc_pipe(i0)%twater_curve
            dfact   =tcurves(twater_curve)%dfact

            do i1=1,nline_g_w
                npairs_wc=wc_pipe(i0)%line_g_w(i1)%npairs_wc
                allocate(value(npairs_wc),heat_wc(npairs_wc))
                pairnode_wc=>wc_pipe(i0)%line_g_w(i1)%pairnode_wc
                value=0.;heat_wc=0.

                if(iwc==1)then
                    ipoin=pairnode_wc(1)
                    jtotv=nodfn(lmdofn(10),ipoin)
                    result_zero(jtotv)=dfact
                    result_first(jtotv)=0.
                elseif(iwc==-1)then
                    ipoin=pairnode_wc(npairs_wc)
                    jtotv=nodfn(lmdofn(10),ipoin)
                    result_zero(jtotv)=dfact
                    result_first(jtotv)=0.
                endif
                do ipairs=1,npairs_wc
                    itotv=nodfn(lmdofn(10),pairnode_wc(ipairs))
                    if(itotv/=0) &
                        value(ipairs)=result_zero(itotv)
                end do

                Ks=>wc_pipe(i0)%line_g_w(i1)%kmatrix_w
                heat_wc=Ks.x.value
                heat_wc=coef1*heat_wc
                !write(7,*)'heat_wc='
                !write(7,30)heat_wc

                do ipairs=1,npairs_wc
                    jtotv=nodfn(lmdofn(10),pairnode_wc(ipairs))
                    tofor(jtotv)=-heat_wc(ipairs)
                end do
                deallocate(value,heat_wc)
                nullify(Ks,pairnode_wc)
            end do
        end do

    endif  !20210417



    if (mdiv==1)then
        refor=refor+(tofor-stfor)
    else
        refor=refor+(toform-stfor)
    endif

    if (mdiv==1)then
        toforx=tofor
    else
        toforx=toform
    endif

    !!int2000
    ! write(7,*)'refor(5218)=',refor(5218)
    !write(7,*)'tofor(17570)=',tofor(17570),'stfor(17570)=',stfor(17570)
    !write(7,*)'refor(17570)=',refor(17570)
    !


    tofory=toforx
    refory=refor
    do itotv=1,ntotv
        nintf=trans(itotv)%nintf

        if (nintf/=0) then
            !write(7,*)'itotv=',itotv,'nintf=',nintf

            toforx(itotv)=0.
            refor(itotv)=0.
            do iintf=1,nintf
                !if (iffix(trans(itotv)%listf(iintf))==0)then  !20211214
                toforx(trans(itotv)%listf(iintf))=  &
                    toforx(trans(itotv)%listf(iintf))+  &
                    tofory(itotv)*trans(itotv)%rintf(iintf)
                refor(trans(itotv)%listf(iintf))=  &
                    refor(trans(itotv)%listf(iintf))+  &
                    refory(itotv)*trans(itotv)%rintf(iintf)
                !if(trans(itotv)%listf(iintf)==5218)write(7,*)trans(itotv)%listf(iintf),itotv,refor(trans(itotv)%listf(iintf))
                !endif  !20211214
            end do
        endif
    end do
    !!int2000
    ! write(7,*)'total force on spring nodes of dam_foundation interface'
    !  do idimn=1,mdofn
    !tfi=0.
    !if(lmdofn(idimn)==0)cycle  !20220302
    !do ipoin=1,npoin
    !    itotv=nodfn(lmdofn(idimn),ipoin) !20220302
    !    if (itotv==0) cycle
    !     if (iffix(itotv)==6) &
    !    tfi=tfi-refor(itotv)
    !end do
    !write(7,*)'idimn=',idimn,'tfi=',tfi
    !  end do

    if(nbackf/=0)then
        ! record the displacements of spring points
        kkdimn=ndimn  !2015/11/17
        if(block_stab==1)kkdimn=3*(ndimn-1) !2015/11/17
        do igapbf=1,nbackf
            igapb=backf(igapbf)%groupb
            npgblock=gapb(igapb)%npgblock
            do i0=1,npgblock
                igaps=gapb(igapb)%nodegblock_igaps(i0)
                ipair=gapb(igapb)%nodegblock_ipairs(i0)
                ij=gapb(igapb)%nodegblock_onetwo(i0)
                ipoin=gaps(igaps)%pairnode(ij,ipair)
                gapb(igapb)%force_ct(1:kkdimn,i0)=-refor(nodfn(1:kkdimn,ipoin))  !坝和地基交界点处位移
            end do
        end do
    endif

    do itotv=1,ntotv  !20191006
        if (iffix(itotv)>0) refor(itotv)=0.
    end do




    if (nflow/=0)then
        do ii=1,nfreeflownode
            ipoin=listfreeflownode(ii)
            itotv=nodfn(lmdofn(8),ipoin)
            if (totveq(itotv)/=0.and.(result_zero(itotv)>.1  &
                .or.(result_zero(itotv)>0..and.flowrate(ipoin)>0.))) then
                refor(itotv)=0.
            endif
        end do
    endif

    do ilink=1,nlinks
        npairs=links(ilink)%npairs
        pairnode=>links(ilink)%pairnode
        link_freedom=>links(ilink)%link_freedom
        do order_freedom=1,cdofn
            if (link_freedom(order_freedom)==1) then
                do ipair=1,npairs
                    first_node   =pairnode(1,ipair)
                    second_node  =pairnode(2,ipair)
                    itotv        =nodfn(order_freedom,first_node)
                    jtotv        =nodfn(order_freedom,second_node)
                    toforx(itotv)  =toforx(itotv)+toforx(jtotv)
                    toforx(jtotv)  =0.0
                    refor (itotv)  =refor (itotv)+refor (jtotv)
                    refor (jtotv)  =0.0
                end do
            endif
        end do
        nullify(pairnode)
    end do

    !write(outact,*)'istep=',istep,'iiter=',iiter,'tofor_ncomn='
    !    allocate(value(1:ndimn))
    !    i0=0
    !       do ipoin=1,npoin
    !          value=0.
    !          do idofn=1,ndimn
    !             itotv=nodfn(idofn,ipoin)
    !             if (itotv/=0)value(idofn)=tofor(itotv)
    !          end do
    !          if(any(abs(value)>.01)) then
    !              i0=i0+1
    !          write(outact,11)i0,ipoin,value
    !          endif
    !       end do
    !       deallocate(value)
11  format(2i10, 6e15.5)
    !12 format(3i10, 3e20.5)

    !write(7,*)'istep=',istep,'iiter=',iiter
    !write(7,*)'ipoin,idofn,itotv,iffix(itotv),toforx(itotv),stfor(itotv),refor(itotv)='
    !
    !       do ipoin=1,npoin
    !          do idofn=1,mdofn  !ndimn
    !           if(lmdofn(idofn)==0)cycle
    !             itotv=nodfn(lmdofn(idofn),ipoin)     !20200220
    !             if(itotv==0)cycle
    !                !if(abs(refor(itotv))>1.e-3) &
    !          write(7,1234)ipoin,idofn,itotv,iffix(itotv),toforx(itotv),stfor(itotv),refor(itotv)
    !           end do
    !       end do
    !
    !1234 format(4I10,3e15.5)
    !
    resid=dot_product(refor,refor)
    retot=dot_product(toforx,toforx)

    max_refor=maxval(abs(refor))
    max_tofor=maxval(abs(toforx))
    print *, 'resid=',resid,'retot=',retot
    resid=sqrt(resid)
    retot=sqrt(retot)
    !   print *, 'resid=',resid,'retot=',retot
    !   print *, 'max_refor=',max_refor,'max_tofor=',max_tofor
    ratio1=resid/retot
    ratio2=max_refor/max_tofor

    write(chkunit,*)'iblks=',iblks,'iincs=',iincs,'istep=',istep,'iiter=',iiter
    write(chkunit,*)'resid=',resid,'retot=',retot
    write(chkunit,10) ratio1
    write(chkunit,*)'max_refor=',max_refor,'max_tofor=',max_tofor
    write(chkunit,20) ratio2
    print *,'ratio1=',ratio1,'ratio2=',ratio2
    !   if(ratio1.ge.1..or.ratio2.ge.1.) stop 'error conver_load'

    deallocate(refor,toforx,refory,tofory)

    print *,'nchek0=',nchek,'toler_force=',toler_force
    if(ratio1>toler_force.or.ratio2>toler_force) nchek=1
    print *,'nchek=',nchek

10  format('ratio for residu norm=',e11.4)
20  format('ratio for residu maxm=',e11.4)
30  format(10e15.5)

    end subroutine conver_load

    subroutine conver_load_w !freq2006

    complex(irk),allocatable:: refor(:)
    real   (irk) resid,retot,max_refor,max_tofor,ratio1,ratio2

    !      ------  Checks convergence for load part

    if (iiter==2)then
        write(7,*)'ndofn='
        do ipoin=1,npoin
            write(7,*)ipoin,nodfn(lmdofn(1:ndimn),ipoin)
        end do
    endif
    print *,'                    istep=',istep,'iiter=',iiter
    print *,'                    conver_load check'
    write(chkunit,*)
    write(chkunit,*)'        Convergence check for Residual force,ttime=,',ttime
    allocate(refor(ntotv))
    nchek=0
    retot=0.0
    refor=(0.0,0.)
    refor=(toforw-stforw)

    do itotv=1,ntotv
        if(iffix(itotv)/=0)then
            refor(itotv)=(0.,0.)
            toforw(itotv)=stforw(itotv)
        endif
    end do

    resid=dot_product(refor,refor)
    retot=dot_product(toforw,toforw)

    max_refor=maxval(abs(refor))
    max_tofor=maxval(abs(toforw))
    print *, 'resid=',resid,'retot=',retot
    resid=sqrt(resid)
    retot=sqrt(retot)
    !   print *, 'resid=',resid,'retot=',retot
    !   print *, 'max_refor=',max_refor,'max_tofor=',max_tofor
    ratio1=resid/retot
    ratio2=max_refor/max_tofor
    !write(chkunit,*)'iiter=',iiter,'stfor,tofor,refor='
    !do itotv=1,ntotv
    !write(chkunit,30)itotv,stforw(itotv),toforw(itotv),refor(itotv)
    !end do


    write(chkunit,*)'iblks=',iblks,'iincs=',iincs,'istep=',istep,'iiter=',iiter
    write(chkunit,*)'resid=',resid,'retot=',retot
    write(chkunit,10) ratio1
    write(chkunit,*)'max_refor=',max_refor,'max_tofor=',max_tofor
    write(chkunit,20) ratio2
    print *,'ratio1=',ratio1,'ratio2=',ratio2
    !   if(ratio1.ge.1..or.ratio2.ge.1.) stop 'error conver_load'

    deallocate(refor)

    if(ratio1>toler_force.or.ratio2>toler_force) nchek=1

10  format('ratio for residu norm=',e11.4)
20  format('ratio for residu maxm=',e11.4)
30  format(i10,10e15.4)

    end subroutine conver_load_w

    subroutine conver_nodal_value

    integer(ink) checki,mdofn1,mdofn2,idimn,itotv,ipoin
    real   (irk),allocatable::tavalue(:),itvalue(:)

    real   (irk) tavnorm,itvnorm,max_tav,max_itv,ratio1,ratio2

    !      ------  Checks convergence for load part

    print *, '                   conver_nodal_value check'
    write(chkunit,*)
    nchek=0

    !!  check convergence for displacement term

    checki=1
    mdofn1=1
    mdofn2=ndimn

    allocate(tavalue(ntotv),itvalue(ntotv))


100 if (any(lmdofn(mdofn1:mdofn2)/=0)) then


        if(checki==1)write(chkunit,*)'              Convergence check for displacement'
        if(checki==2)write(chkunit,*)'              Convergence check for rotations'
        if(checki==3)write(chkunit,*)'              Convergence check for pressure(P)'
        if(checki==4)write(chkunit,*)'              Convergence check for pressure(Pw)'
        if(checki==5)write(chkunit,*)'              Convergence check for pressure(Pa)'
        if(checki==6)write(chkunit,*)'              Convergence check for temperature'
        if(checki==1)print *,'              Convergence check for displacement'
        if(checki==2)print *,'              Convergence check for rotations'
        if(checki==3)print *,'              Convergence check for pressure(P)'
        if(checki==4)print *,'              Convergence check for pressure(Pw)'
        if(checki==5)print *,'              Convergence check for pressure(Pa)'
        if(checki==6)print *,'              Convergence check for temperature'
        tavalue=0.0
        itvalue=0.0
        do idimn=mdofn1,mdofn2
            if (lmdofn(idimn)/=0) then
                do ipoin=1,npoin
                    itotv=nodfn(lmdofn(idimn),ipoin)
                    if (itotv>0)then
                        tavalue(itotv)=result_zero(itotv)
                        itvalue(itotv)=delitfi(itotv)
                        !if(iiter>1.and.abs(itvalue(itotv)/tavalue(itotv))>1.e-3)then
                        !write(7,*)'itotv=',itotv,'itvalue=',itvalue(itotv),'tavalue=',tavalue(itotv)
                        !endif
                    endif
                end do
            endif
        end do

        !write(7,*)'tavalue,itvalue'
        !do itotv=1,ntotv
        !    write(7,*)itotv,tavalue(itotv),itvalue(itotv)
        !end do

        tavnorm=dot_product(tavalue,tavalue)
        itvnorm=dot_product(itvalue,itvalue)
        max_tav=maxval(abs(tavalue))
        max_itv=maxval(abs(itvalue))

        tavnorm=sqrt(tavnorm)
        itvnorm=sqrt(itvnorm)
        if (max_tav.le.1.e-10.or.tavnorm.le.1.e-10)then
            nchek=0
            goto 101
        endif
        !      print *,'tatfn=',tavnorm,'itifn=',itvnorm
        write(7,*)'max_tatf=',max_tav,'max_itif=',max_itv
        ratio1=itvnorm/tavnorm
        ratio2=max_itv/max_tav

        write(chkunit,*)'iiter=',iiter,'tatfn=',tavnorm,'itifn=',itvnorm
        write(chkunit,10) ratio1
        write(chkunit,*)'max_tatf=',max_tav,'max_itif=',max_itv
        write(chkunit,20) ratio2
        print *,'ratio1=',ratio1,'ratio2=',ratio2

        if(ratio1>toler_var(mdofn1).or.ratio2>toler_var(mdofn1)) nchek=1

    end if   !! for convergence checki


    if (nchek==1) then
        print *,'                           not converged for checki=',checki
        write(chkunit,*)'                   not converged for checki=',checki
        deallocate(tavalue,itvalue)
        return
    endif

101 continue


    if (checki==6) then
        deallocate(tavalue,itvalue)
        return

    else

        checki=checki+1
        if(checki==2)             mdofn1=4
        if(checki==2.and.ndimn==2)mdofn2=4
        if(checki==2.and.ndimn==3)mdofn2=6

        if(checki==3)mdofn1=7
        if(checki==3)mdofn2=7

        if(checki==4)mdofn1=8
        if(checki==4)mdofn2=8

        if(checki==5)mdofn1=9
        if(checki==5)mdofn2=9

        if(checki==6)mdofn1=10
        if(checki==6)mdofn2=10

        if (mdofn<mdofn1.or.(checki==5.and.(outintr.ne.0.or.outinp/=0)).or.(checki==3.and.outinp==8))then
            print *, 'nchek=',nchek
            deallocate(tavalue,itvalue)
            return
        endif
        goto 100
    endif


10  format('ratio for norm=',e11.4)
20  format('ratio for maxm=',e11.4)

    end subroutine conver_nodal_value
