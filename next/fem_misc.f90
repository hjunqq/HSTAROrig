    subroutine write_stiff_u

    integer(ink) igroup,ielgroup,ielem,matno,icreep
    character(10)field1,name
    character(20)material

    rewind(midstif)
    DO igroup =1,ngroup
        field1= group(igroup)%fieldid
        if (field1(1:1)=='U'.and.appear(igroup)>0)  then
            matno = group(igroup)%matno
            name=props(matno)%name
            material=props(matno)%mechanical%solid%material
            icreep =props(matno)%mechanical%solid%icreep
            if (name=='CONTACT'.or.material/='ELASTIC_ISOTROPIC'.or.   &
                (material=='ELASTIC_ISOTROPIC'.and.icreep/=0))then
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    write(midstif) element(ielem)%field(1)%khandmc(1)%fstif
                end do
            endif
        endif
    end do

    end subroutine write_stiff_u

    SUBROUTINE neuman_expan

    !! need the zero, firt, second deritives of the result, store in the ntotv order
    !!          result_zero, result_first, result_second

    character(10)fieldid,name
    character(20)material
    integer(ink) igroup,index,nevab,ielgroup,ielem,itotv,iter_neuman,ichk,eigen,matno,icreep
    integer(ink), pointer::ldofs(:)
    real   (irk), allocatable::value(:),fstif0(:,:),refu(:),result0(:),resultn(:),deliu(:), &
        eload(:),result_org(:)
    real   (irk), pointer::fstif(:,:)
    real   (irk) err,fmfold,rld
    allocate(refu(ntotv),result0(ntotv),resultn(ntotv),deliu(ntotv))
    allocate(result_org(ntotv))
    eigen=0
    result_org=result
10  iter_neuman=0
    fmfold=1.
    result0=result/fmfold
    resultn=result/fmfold
1   iter_neuman=iter_neuman+1
    where(iffix==0)
        deliu=result0
    elsewhere
        deliu=0.
    endwhere
    refu=0.0
    rewind(midstif)
    DO igroup =1,ngroup
        fieldid=group(igroup)%fieldid
        if (appear(igroup)>0.and.fieldid(1:1)=='U') then
            matno=group(igroup)%matno
            name=props(matno)%name
            material=props(matno)%mechanical%solid%material
            icreep =props(matno)%mechanical%solid%icreep
            if (name=='CONTACT'.or.material/='ELASTIC_ISOTROPIC'.or.   &
                (material=='ELASTIC_ISOTROPIC'.and.icreep/=0))then
                ! get information from the group level
                index=group(igroup)%index
                nevab=elkn(index)%el_field(1)%nnode_f*group(igroup)%dof(1)%nfdof
                allocate(value(nevab),fstif0(nevab,nevab),eload(nevab))
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    read(midstif)fstif0
                    fstif=>element(ielem)%field(1)%khandmc(1)%fstif
                    ldofs=>element(ielem)%field(1)%ldofs_f
                    fstif0=fstif/fmfold-fstif0
                    value=deliu(ldofs)
                    eload=MATMUL(fstif0,value)
                    refu(ldofs)=refu(ldofs)+eload
                    nullify(fstif,ldofs)
                end do       !!ielgroup
                deallocate(value,fstif0,eload)
            endif
        endif
    enddo !for igroup

    rvector=0.
    do itotv=1,ntotv
        if(totveq(itotv)/=0)     &
            rvector(totveq(itotv))=rvector(totveq(itotv))+refu(itotv)
    end do

    operation='SOLVE'
    call solve
    resultn=resultn+(-1)**iter_neuman*result
    err=MAXVAL(abs(result))/MAXVAL(abs(resultn))
    ichk=0
    write(chkunit,*)'iter_neuman==',iter_neuman,'err=',err
    if(err.le.1.e-1.or.iter_neuman.eq.30)ichk=1
    if (ichk==0)then
        !      if((iter_neuman==20.or.err.ge.1.0).and.eigen==0) then
        if (err.ge.1.0.and.eigen==0) then
            print *,'neuman expansion failed!'
            stop
            eigen=1
            call eigv(rld)
            fmfold=(rld+1.)/2.+.2
            result=result_org
            goto 10
        else
            result0=result
            goto 1
        endif
    endif
    print *,'iter_neuman=',iter_neuman,'err=',err
    result=resultn
    deallocate(refu,result0,resultn,deliu,result_org)

    END SUBROUTINE neuman_expan

    SUBROUTINE kdelt(ic,ics,result0,refu)

    character(10)fieldid
    integer(ink) igroup,index,ic,ics,nevab, ielgroup, ielem, itotv
    integer(ink),pointer::ldofs(:)
    real   (irk),allocatable::value(:),fstif0(:,:), eload(:),deliu(:)
    real   (irk),pointer::fstif(:,:)
    real   (irk) result0(:),refu(:)

    allocate(deliu(ntotv))
    where(iffix==0)
        deliu=result0
    elsewhere
        deliu=0.
    endwhere
    refu=0.0
    if(ic==2)rewind(midstif)
    DO igroup =1,ngroup
        fieldid = group(igroup)%fieldid
        if (appear(igroup)>0.and.fieldid(1:1)=='U') then
            ! get information from the group level
            index=group(igroup)%index
            nevab=elkn(index)%el_field(1)%nnode_f*group(igroup)%dof(1)%nfdof
            allocate(value(nevab),fstif0(nevab,nevab),eload(nevab))
            DO ielgroup=1,group(igroup)%nelgroup
                ielem=group(igroup)%list(ielgroup)
                if(ic==2)read(midstif)fstif0
                fstif=>element(ielem)%field(1)%khandmc(1)%fstif
                ldofs=>element(ielem)%field(1)%ldofs_f
                if (ic==1)then
                    fstif0=fstif
                else
                    fstif0=fstif-fstif0
                endif
                value=deliu(ldofs)
                eload=MATMUL(fstif0,value)
                refu(ldofs)=refu(ldofs)+eload
                nullify(fstif,ldofs)
            enddo !!ielgroup
            deallocate(value,fstif0,eload)
        endif
    enddo     !!  for igroup
    deallocate(deliu)
    if(ics==0)return
    rvector=0.
    do itotv=1,ntotv
        if(totveq(itotv)/=0)rvector(totveq(itotv))=rvector(totveq(itotv))+refu(itotv)
    end do
    operation='SOLVE'
    call solve

    END SUBROUTINE kdelt

    SUBROUTINE EIGV (rld)

    integer(ink) num
    real   (irk) rec,aa,bb,rld,re,r3,r4
    real   (irk), allocatable::refu(:)
    print *,'******','EIGENNALUE CALCULATION','*******'
    allocate(refu(ntotv))
    NUM=0
    BB=0.0
    REC=1.0E-7
    result=1.
12  call kdelt(2,1,result,refu)
    RLD=maxval(abs(result))
    IF(rld.LT.0.00001) print *, 'AAA1'
    result=result/rld
    RE=0.0
    R3=0.0
    R4=0.0
    call kdelt(1,0,result,refu)
    r3=result.d.refu
    call kdelt(2,0,result,refu)
    r4=result.d.refu
    RLD=abs(R4/R3)
    AA=RLD
    RE=ABS(AA-BB)/AA
    IF (RE.GT.REC) BB=AA
    NUM=NUM+1
    IF (NUM.GT.150) GOTO 165
    print *, 'num=',num,'rld=',rld
    IF (RE.LE.REC) GOTO 165
    GO TO 12
165 CONTINUE
    !   WO=1.0/SQRT(AA)
    !   TO=6.2832/WO
    !   FO=1./TO
    deallocate(refu)
    WRITE(chkunit,*) 'num=',num,'rld=',rld

    RETURN

    END subroutine eigv

    SUBROUTINE RESTA_READ_WRITE(ic)

    character(10) field1,name
    integer(ink) igroup,index,order_int,ielem,matno,ic,ielgroup,icreep,igapb,nrdof,igaps
    integer(ink) kinit_g !20211214
    rewind(restaunit)
    if (ic==-1) then
        write(restaunit)iblks,iincs,ttime,line_load_block(1:iblks),line_temp_block(1:iblks),trstep
        write(restaunit)group(1:ngroup)%btime,group(1:ngroup)%ditime_1
        if(allocated(result_zero))  write(restaunit)result_zero
        if(allocated(result_first)) write(restaunit)result_first
        if(allocated(result_second))write(restaunit)result_second
        if(allocated(prstat))       write(restaunit)prstat
        if(allocated(torel))        write(restaunit)torel
        if(allocated(toforl))       write(restaunit)toforl
        if(allocated(fexta))        write(restaunit)fexta  !!nstoks

        DO igroup =1,ngroup
            field1= group(igroup)%fieldid
            if (field1(1:1)=='U')  then
                kinit_g=group(igroup)%kinit_g
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    write(restaunit)element(ielem)%field(1)%gpvar
                    write(restaunit)element(ielem)%field(1)%tload
                    write(restaunit)element(ielem)%field(1)%eload
                    if(kinit_g==2)write(restaunit)element(ielem)%stres0
                end do
            endif
        end do
        DO igroup =1,ngroup
            field1= group(igroup)%fieldid
            index = group(igroup)%index
            if (field1(1:1)=='U')  then
                matno = group(igroup)%matno
                name  = props(matno)%name
                icreep=props(matno)%mechanical%solid%icreep
                !! contact
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    if (name=='CONTACT') then
                        write(restaunit)element(ielem)%field(1)%gapg
                        write(restaunit)element(ielem)%field(1)%gapn
                        write(restaunit)element(ielem)%field(1)%state
                    endif
                    if (icreep==2) then
                        write(restaunit)element(ielem)%field(1)%omega
                        write(restaunit)element(ielem)%field(1)%dsig
                    end if
                end do
            endif
            !!contact
            if (field1(1:2)=='UW') then
                if (name=='NSSoilPZ') then
                    order_int=elkn(index)%el_field(1)%order_intrules(1)
                    DO ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        write(restaunit)element(ielem)%egaus(order_int)%iload
                        write(restaunit)element(ielem)%egaus(order_int)%vdval
                    end do
                endif
            endif
        end do

        do igaps=1,ngaps
            write(restaunit)gaps(igaps)%state
            write(restaunit)gaps(igaps)%ctforce
            if(kinit==2)write(restaunit)gaps(igaps)%ctforce_stres0  !20220623
            write(restaunit)gaps(igaps)%gap
            write(restaunit)gaps(igaps)%dxyz
        end do

        do igapb=1,ngapb
            nrdof=gapb(igapb)%nrdof
            if(nrdof==0)cycle
            write(restaunit)gapb(igapb)%rdisp_zero
            if(type_problem=='F')then
                write(restaunit)gapb(igapb)%rdisp_first
                write(restaunit)gapb(igapb)%rdisp_second
            endif
        end do


        print *,'***********RESTART FILE IS UPDATED!*******'
    else if(ic==1) then
        Print *,'restart reading file'
        read(restaunit)lblks,lincs,lttime,line_load_block(1:lblks),line_temp_block(1:lblks),trstep
        print *,'lblks=',lblks,'lincs=',lincs
        read(restaunit)group(1:ngroup)%btime,group(1:ngroup)%ditime_1
        if(allocated(result_zero))  read(restaunit)result_zero
        if(allocated(result_first)) read(restaunit)result_first
        if(allocated(result_second))read(restaunit)result_second
        if(allocated(prstat))       read(restaunit)prstat
        if(allocated(torel))        read(restaunit)torel
        if(allocated(toforl))       read(restaunit)toforl
        if(allocated(fexta))        read(restaunit)fexta  !!nstoks
        DO igroup =1,ngroup
            field1= group(igroup)%fieldid
            if (field1(1:1)=='U')  then
                kinit_g=group(igroup)%kinit_g  !20211214
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    read(restaunit)element(ielem)%field(1)%gpvar0
                    read(restaunit)element(ielem)%field(1)%tload
                    read(restaunit)element(ielem)%field(1)%eload
                    element(ielem)%field(1)%gpvar=element(ielem)%field(1)%gpvar0
                    if(kinit_g==2)read(restaunit)element(ielem)%stres0
                end do
            endif
        end do
        DO igroup =1,ngroup

            field1= group(igroup)%fieldid
            index = group(igroup)%index
            !! contact
            if (field1(1:1)=='U')  then
                matno = group(igroup)%matno
                name=props(matno)%name
                icreep=props(matno)%mechanical%solid%icreep
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    if (name=='CONTACT') then
                        read(restaunit)element(ielem)%field(1)%gapg0
                        read(restaunit)element(ielem)%field(1)%gapn0
                        read(restaunit)element(ielem)%field(1)%state0
                        element(ielem)%field(1)%gapg=element(ielem)%field(1)%gapg0
                        element(ielem)%field(1)%gapn=element(ielem)%field(1)%gapn0
                    endif
                    if (icreep==2) then
                        read(restaunit)element(ielem)%field(1)%omega
                        read(restaunit)element(ielem)%field(1)%dsig
                    end if
                end do
            endif
            !!contact
            if (field1(1:2)=='UW') then
                if (name=='NSSoilPZ') then
                    order_int=elkn(index)%el_field(1)%order_intrules(1)
                    DO ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        read(restaunit)element(ielem)%egaus(order_int)%iload0
                        read(restaunit)element(ielem)%egaus(order_int)%vdval0
                    end do
                endif
            endif
        end do


        do igaps=1,ngaps
            read(restaunit)gaps(igaps)%state
            read(restaunit)gaps(igaps)%ctforce
            if(kinit==2)read(restaunit)gaps(igaps)%ctforce_stres0  !20220623
            read(restaunit)gaps(igaps)%gap
            read(restaunit)gaps(igaps)%dxyz

            gaps(igaps)%state0=gaps(igaps)%state
            gaps(igaps)%ctforce0=gaps(igaps)%ctforce
            gaps(igaps)%gap0=gaps(igaps)%gap
            gaps(igaps)%dxyz0=gaps(igaps)%dxyz
        end do

        do igapb=1,ngapb
            nrdof=gapb(igapb)%nrdof
            if(nrdof==0)cycle
            read(restaunit)gapb(igapb)%rdisp_zero
            if(type_problem=='F')then
                read(restaunit)gapb(igapb)%rdisp_first
                read(restaunit)gapb(igapb)%rdisp_second
            endif
        end do

        print *,'***********RESTART FILE IS READ!*******'
    endif


    END SUBROUTINE RESTA_READ_WRITE

    subroutine safety_factor
    character(20)material,criteria,name
    integer(ink) ielem,igroup,iforce,jgroup,matno,ielgroup
    integer(ink) index,ngaus,igaus,order_int,nstre,npairs,igaps,ipairs,isafety,istate
    integer(ink),allocatable::iii(:)
    real   (irk),allocatable::ftang(:),fresi(:),ftang_gaps(:),fresi_gaps(:),sfactor(:),ps(:)  !20200411
    real   (irk) ft,fn,uniax,frict,k_safety,elcod_local,yld,aera, dilan,sigma0,ftangt,fresit,t1,t2,tt,areat,damage

    if(nforce==0.and.ngaps==0) return
    if(nforce==0) goto 10
    allocate(ftang(nforce),fresi(nforce),iii(ngroup))
    allocate(ps(ndimn))
    iii=0
    ftang=0.
    fresi=0.
    areat = 0.

    do iforce=1,nforce
        if(nforce_appear(iforce)/=1.and.nforce_appear(iforce)/=3)cycle !zhao 2010
        do jgroup=1,surface_force(iforce)%lgroup
            igroup=surface_force(iforce)%list(jgroup)
            iii(igroup)=1
            index = group(igroup)%index
            nstre =group(igroup)%nstre
            matno =group(igroup)%matno
            name=props(matno)%name
            select case(trim(props(matno)%mechanical%solid%material))
            case('CLASSICALEP')
                order_int=elkn(index)%el_field(1)%order_intrules(1)
                ngaus = elkn(index)%ggaus(order_int)%ngaus
                elcod_local=group(igroup)%elcod_local
                uniax   =props(matno)%mechanical%solid%classicalEP%sigma0
                frict   =props(matno)%mechanical%solid%classicalEP%frict_angle
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    do igaus=1,ngaus
                        yld=element(ielem)%field(1)%gpvar(nstre+2,igaus)
                        if (yld/=2.) then
                            aera=1./elcod_local*element(ielem)%egaus(1)%djacb(igaus)
                            areat = areat + aera
                            fn=element(ielem)%field(1)%ntstress(1,igaus)
                            ft=element(ielem)%field(1)%ntstress(2,igaus)
                            fresi(iforce)=fresi(iforce)-fn*aera*tand(frict)+uniax*aera
                            ftang(iforce)=ftang(iforce)+ft*aera
                        endif
                    end do
                end do

            case('GOODMAN')

                order_int=elkn(1)%el_field(1)%order_intrules(1)
                ngaus=elkn(1)%ggaus(order_int)%ngaus
                elcod_local=group(igroup)%elcod_local
                frict = props(matno)%mechanical%solid%Goodman%frict_angle
                uniax = props(matno)%mechanical%solid%Goodman%uniax_cohes
                if (name=="CONTACT")then
                    do ielgroup = 1, group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        do igaus = 1, ngaus
                            ps = element(ielem)%field(1)%gpvar(1:ndimn,igaus)
                            ps = element(ielem)%field(1)%gpvar(1:ndimn,igaus)
                            fn = ps(2)
                            ft = ps(1)
                            if (element(ielem)%field(1)%state(igaus)=='contact') then
                                aera=element(ielem)%aera_local(igaus)
                                areat = areat + aera
                                fresi(iforce)=fresi(iforce)-fn*aera*tand(frict)+uniax*aera
                            else
                                fn = ps(2)
                                fresi(iforce)=fresi(iforce)-fn*aera*tand(frict)
                            endif
                            ftang(iforce)=ftang(iforce)+ft*aera
                        enddo
                    enddo
                else
                    do ielgroup = 1, group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        do igaus =1, ngaus
                            ps = element(ielem)%field(1)%gpvar(1:ndimn,igaus)
                            fn = ps(2)
                            ft = ps(1)
                            if (ps(2)<1e-10) then
                                aera=element(ielem)%aera_local(igaus)
                                areat = areat + aera
                                fresi(iforce)=fresi(iforce)-fn*aera*tand(frict)+uniax*aera
                            else
                                fn = ps(2)
                                fresi(iforce)=fresi(iforce)-fn*aera*tand(frict)
                            endif
                            ftang(iforce)=ftang(iforce)+ft*aera
                        enddo
                    enddo
                endif

            end select

        end do
        write(chkunit,*)'iforce=',iforce
        write(chkunit,*)'fresi=',fresi(iforce),'ftang=',ftang(iforce)
    end do !iforce

    fresit=0.;ftangt=0.
    fresit=sum(fresi)
    ftangt=sum(ftang)
    K_safety=kstab*fresit/ftangt
    !  kstab=K_safety
    write(chkunit,*)'k_safety=',k_safety, 'area_sum=',areat

10  continue
    if(ngaps/=0)then
        allocate(fresi_gaps(ngaps),ftang_gaps(ngaps))
        if(nsafety_gaps>=0)allocate(sfactor(nsafety_gaps))
        do isafety=1,nsafety_gaps  !isafety  20200409
            fresi_gaps=0.;ftang_gaps=0.
            do igaps=1,ngaps
                if(safety_gaps_appear(igaps,isafety)==0)cycle
                npairs=gaps(igaps)%npairs
                !write(7,*)'igaps=',igaps,'fxyz=',sum(gaps(igaps)%ctforce(1,:)),sum(gaps(igaps)%ctforce(2,:))
                do ipairs=1,npairs
                    if(gaps(igaps)%state(ipairs)==0)cycle  !20211005

                    !write(7,*)'ipairs=',ipairs,'gaps(igaps)%cohes(ipairs)=',gaps(igaps)%cohes(ipairs),'gaps(igaps)%frict(ipairs) =',gaps(igaps)%frict(ipairs)
                    fresi_gaps(igaps)=fresi_gaps(igaps)+gaps(igaps)%aera(ipairs)*gaps(igaps)%cohes(ipairs)-gaps(igaps)%ctforce(ndimn,ipairs)*gaps(igaps)%frict(ipairs)
                    t1=gaps(igaps)%ctforce(1,ipairs)
                    if (ndimn==3)t2=gaps(igaps)%ctforce(2,ipairs)
                    tt=abs(t1)
                    if (ndimn==3)tt=sqrt(t1**2+t2**2)
                    ftang_gaps(igaps)=ftang_gaps(igaps)+tt
                    !write(7,*)'sfl=',(gaps(igaps)%aera(ipairs)*gaps(igaps)%cohes(ipairs)-gaps(igaps)%ctforce(ndimn,ipairs)*gaps(igaps)%frict(ipairs))/tt
                end do
            enddo
            fresit=0.;ftangt=0
            fresit=fresit+sum(fresi_gaps)
            ftangt=ftangt+sum(ftang_gaps)
            write(chkunit,*)'isafety=',isafety
            write(chkunit,*)'fresit=',fresit,'ftangt=',ftangt
            ! write(chkunit,*)'fresi_gaps=',fresi_gaps
            !write(chkunit,*)'ftang_gaps=',ftang_gaps

            K_safety=kstab*fresit/ftangt
            sfactor(isafety)=K_safety
            write(chkunit,*)'k_safety=',k_safety
        end do  !isafety  20200409
        write(7,*)'minimum safety=',minval(sfactor)
        deallocate(sfactor)
    end if



    if(nforce/=0) &
        deallocate(ftang,fresi)
    if(ngaps/=0) &
        deallocate(ftang_gaps,fresi_gaps)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! zhao 2006/03/09
    if(nforce/=0) then
        do igroup=1,ngroup
            if(iii(igroup)/=1)cycle
            matno =group(igroup)%matno
            material=props(matno)%mechanical%solid%material
            if(material/='CLASSICALEP')cycle
            criteria=props(matno)%mechanical%solid%ClassicalEP%criteria
            if (kstab/=0.and.criteria(1:2)=='MC') then
                frict=props(matno)%mechanical%solid%ClassicalEP%frict_angle_ini
                dilan=props(matno)%mechanical%solid%ClassicalEP%dilan_angle_ini
                frict=tand(frict)/kstab
                dilan=tand(dilan)/kstab
                frict=atand(frict)
                dilan=atand(dilan)
                sigma0=props(matno)%mechanical%solid%ClassicalEP%sigma0_ini
                sigma0=sigma0/kstab
                props(matno)%mechanical%solid%ClassicalEP%frict_angle=frict
                props(matno)%mechanical%solid%ClassicalEP%dilan_angle=dilan
                props(matno)%mechanical%solid%ClassicalEP%sigma0=sigma0
            endif
        enddo
        deallocate(iii)
    endif

    end subroutine safety_factor

    subroutine FORCE_interface
    character(20)material,criteria
    integer(ink) ielem,idofn,inode,idimn,igroup,iforce,iextrf, &
        jelem,jgroup,ipoin,ne_unode,ipface,matno,ig,mgroup,index,nevab
    integer(ink),allocatable::ice(:),icp(:)
    integer(ink),pointer::lnods(:),liste(:)
    real   (irk),pointer::eload(:),tload(:),elcod(:,:)

    allocate(ice(nelem)) !,icp(npoin))
    do iforce=1,nforce
        ice=0
        liste=>surface_force(iforce)%liste
        ice(liste)=1
        if (nliste==2)then
            nullify(liste)
            liste=>surface_force(iforce)%liste1
            ice(liste)=1
        endif
        nullify(liste)
        surface_force(iforce)%ftfor=0.0

        !npface=surface_force(iforce)%npface
        !icp=0
        !do ipface=1,npface
        !   ipoin=surface_force(iforce)%list_npface(ipface)
        !	icp(ipoin)=ipface
        !enddo

        !       do ielem=1,nelem
        !	      if(ice(ielem)==0)cycle
        !            lnods=>element(ielem)%field(1)%lnods_f
        !            eload=>element(ielem)%field(1)%eload
        !            tload=>element(ielem)%field(1)%tload
        !			do inode=1,size(lnods)
        !			   ipoin=lnods(inode)
        !			   if(icp(ipoin)==0)cycle
        !               do idimn=1,ndimn
        !                  idofn=(inode-1)*ndimn+idimn
        !                  surface_force(iforce)%ftfor(idimn,icp(ipoin))=   &
        !                  surface_force(iforce)%ftfor(idimn,icp(ipoin))+eload(idofn)-tload(idofn)
        !                end do
        !			enddo
        !			nullify(eload,tload)
        !	   enddo



        do jgroup=1,surface_force(iforce)%lgroup
            igroup=surface_force(iforce)%list(jgroup)
            index=group(igroup)%index
            if (appear(igroup).gt.0) then
                matno = group(igroup)%matno
                material=props(matno)%mechanical%solid%material
                npface=surface_force(iforce)%npface
                do ipface=1,npface
                    ipoin=surface_force(iforce)%list_npface(ipface)
                    mgroup=listp_group(ipoin)%mgroup
                    do ig=1,mgroup
                        if(listp_group(ipoin)%listg(ig)==igroup) goto 1
                    end do
                    goto 2
1                   ne_unode=group(igroup)%unode(listp_group(ipoin)%listp(ig))%ne_unode

                    do jelem=1,ne_unode
                        ielem=group(igroup)%unode(listp_group(ipoin)%listp(ig))%list(jelem)
                        if (ice(ielem)==1) then

                            elcod=>element(ielem)%field(1)%elcod_f

                            lnods=>element(ielem)%field(1)%lnods_f
                            eload=>element(ielem)%field(1)%eload
                            tload=>element(ielem)%field(1)%tload
                            nevab=size(eload)
                            do inode=1,size(lnods)
                                if (lnods(inode)==ipoin) then
                                    do idimn=1,ndimn
                                        idofn=(inode-1)*ndimn+idimn
                                        if(index==20)idofn=(inode-1)*nevab/2+idimn !steel 2008
                                        surface_force(iforce)%ftfor(idimn,ipface)=   &
                                            surface_force(iforce)%ftfor(idimn,ipface)+eload(idofn)-tload(idofn)
                                    end do
                                end if
                            end do
                            nullify(tload,eload,lnods)
                            nullify(elcod)
                        end if
                    end do !jelem
2                   continue
                end do ! ipface
            endif  !for appear(igroup).gt.0
        end do !jgroup

        !nullify(liste)

        if (nextrf/=0)then  !2004/9/11
            if (istep>nextrf)then
                do iextrf=1,nextrf-1
                    surface_force(iforce)%ftfor_ext(:,:,iextrf+1)=  &
                        surface_force(iforce)%ftfor_ext(:,:,iextrf)
                end do
                surface_force(iforce)%ftfor_ext(:,:,1)=surface_force(iforce)%ftfor
            else
                if(istep==1)surface_force(iforce)%ftfor_ext(:,:,:)=0.
                do iextrf=1,istep-1
                    surface_force(iforce)%ftfor_ext(:,:,iextrf+1)=  &
                        surface_force(iforce)%ftfor_ext(:,:,iextrf)
                end do
                surface_force(iforce)%ftfor_ext(:,:,1)=surface_force(iforce)%ftfor
            endif
        endif   !2004/9/11
    end do !iforce
    deallocate(ice) !,icp)

    end subroutine FORCE_interface

    subroutine write_force_interface
    integer(ink) ipoin,npface,iforce,jpoin
    real   (irk),allocatable::force0(:)
    allocate(force0(ndimn))
    write(ftfunit,*)'ttime=',ttime

    do iforce=1,nforce
        if(nforce_appear(iforce)/=2.and.nforce_appear(iforce)/=3)cycle !zhao 2010
        write(ftfunit,*)'iforce=',iforce
        npface=surface_force(iforce)%npface
        do ipoin=1,npface
            jpoin=surface_force(iforce)%list_npface(ipoin) !steel 2008
            force0=surface_force(iforce)%ftfor(:,ipoin)

            if(icpnorm(jpoin)/=0)surface_force(iforce)%ftfor(:,ipoin)=transpose(prot(:,:,jpoin)).x.force0

            write(ftfunit,1)ipoin,surface_force(iforce)%list_npface(ipoin),surface_force(iforce)%ftfor(:,ipoin)

        end do
        write(ftfunit,*)'total force on x =',sum(surface_force(iforce)%ftfor(1,1:npface))
        write(ftfunit,*)'total force on y =',sum(surface_force(iforce)%ftfor(2,1:npface))
        if(ndimn==3)write(ftfunit,*)'total force on z =',sum(surface_force(iforce)%ftfor(3,1:npface))
    end do
1   format(1x,2i10,3e20.4)

    deallocate(force0)

    end subroutine write_force_interface

    subroutine judge_fine_mesh(igroup,icjr)

    integer(ink) igroup,icjr,nstre,ngaus,igaus,index,ielem,ielgroup
    real   (irk) sx,sy,sxy,delta
    real   (irk),allocatable::smain(:),rr(:,:),stres(:)

    index = group(igroup)%index
    ngaus = elkn(index)%ggaus(1)%ngaus
    nstre=4
    if(ndimn==3)nstre=6
    allocate(stres(nstre),rr(ndimn,ndimn),smain(ndimn))

    icjr=0
    DO ielgroup = 1,group(igroup)%nelgroup
        ielem    = group(igroup)%list(ielgroup)
        do igaus=1,ngaus
            stres=element(ielem)%field(1)%gpvar(1:nstre,igaus)
            if (ndimn==2) then
                sx=stres(1)
                sy=stres(2)
                sxy=stres(3)
                delta=sqrt((sx-sy)**2/4+sxy**2)
                smain=0.
                if(delta.lt.1.e-5) goto 12
                smain(1)=(sx+sy)/2.+delta
                smain(2)=(sx+sy)/2.-delta
12              continue
            else
                call stresmr ( stres, smain, rr)
            endif
            if(smain(1)>500.)icjr=1
        end do
    end do
    deallocate(stres,rr,smain)

    end subroutine judge_fine_mesh

    subroutine find_remesh_element1

    integer(ink) ielem,nstre
    integer(ink),allocatable::rele(:)
    allocate(rele(nelem))
    rele=0
    !if(iblks>1)return
    do ielem=1,nelem
        if(ice0(ielem)==1) goto 10
        nstre=4
        if(ndimn==3)nstre=6
        if(any(element(ielem)%field(1)%gpvar(nstre+3,:)>valv1))rele(ielem)=1
10      continue
    end do
    nelc=sum(rele)
    if (nelc>0)then
        allocate(listnelc(nelc))
        nelc=0
        do ielem=1,nelem
            if (rele(ielem)==1)then
                nelc=nelc+1
                listnelc(nelc)=ielem
            endif
        end do
    endif
    deallocate(rele)

    end subroutine find_remesh_element1

    subroutine find_remesh_element2

    integer(ink) ielem,nstre
    integer(ink),allocatable::rele(:)
    allocate(rele(nelem1))
    rele=0
    !return
    do ielem=1,nelem1
        if(jce1(ielem)==1) goto 10
        nstre=4
        if(ndimn==3)nstre=6
        if(any(element1(ielem)%field(1)%gpvar(nstre+3,:)>valv2))rele(ielem)=1
10      continue
    end do
    nelc1=sum(rele)
    if (nelc1>0)then
        allocate(listnelc1(nelc1))
        nelc1=0
        do ielem=1,nelem1
            if (rele(ielem)==1)then
                nelc1=nelc1+1
                listnelc1(nelc1)=ielem
            endif
        end do
    endif
    deallocate(rele)

    end subroutine find_remesh_element2

    SUBROUTINE  find_remesh_stran

    !*********************************************************************
    !
    !*** by aera weighting average (only for Q4 and B8) elements
    !
    !********************************************************************
    character(10)fieldid,class,material,name
    integer(ink) igroup,index,order_int,ngaus,nevab,idimn,  &
        ipoin,ielem,matno,nnode,ielgroup,nstre
    integer(ink), pointer::lnods(:),ldofs(:)
    integer(ink),allocatable::rele(:)
    real   (irk) elcod_local,stranmax,stranmax1,tstran
    real   (irk),allocatable::stran(:),cartd(:,:),eldis(:),valun(:),valun1(:)

    allocate(valun(nelem))
    if(nelem1>0)allocate(valun1(nelem1))
    valun=0.
    if(nelem1>0)valun1=0.
    do igroup=1,ngroup

        fieldid=group(igroup)%fieldid
        if (appear(igroup)>0.and.fieldid(1:1)=='U')then
            class=group(igroup)%class

            matno = group(igroup)%matno
            name=props(matno)%name
            material=props(matno)%mechanical%solid%material
            index=group(igroup)%index
            nnode=elkn(index)%nnode
            elcod_local=group(igroup)%elcod_local
            if (elcod_local==0..and.fieldid(1:1)=='U'.and.class=='CO'.and.material/='GOODMAN'  &
                .and.name/='CONTACT'.and.nnode/=2)then
                order_int=elkn(index)%el_field(1)%order_intrules(1)

                ngaus=elkn(index)%ggaus(order_int)%ngaus
                nstre=group(igroup)%nstre
                nevab=nnode*ndimn
                allocate (stran(ndimn),cartd(ndimn,nnode),eldis(nevab))
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    if (ice0(ielem)==0)then
                        lnods=>element(ielem)%field(1)%lnods_f
                        ldofs=>element(ielem)%field(1)%ldofs_f
                        eldis = result_zero(ldofs)
                        valun(ielem)=0.
                        !                  do igaus=1,ngaus
                        !                     cartd=element(ielem)%egaus(order_int)%cartd(:,:,igaus)
                        !                     do idimn=1,ndimn
                        !                        stran(idimn)=0.
                        !                        do inode=1,nnode
                        !                           stran(idimn)=stran(idimn)+cartd(idimn,inode)*eldis(ndimn*(inode-1)+idimn)
                        !                        end do
                        !                     end do
                        !                     tstran=sqrt(sum(stran(1:ndimn)**2))
                        !                     valun(ielem)=valun(ielem)+tstran
                        !                  end do
                        valun(ielem)=sum(element(ielem)%field(1)%gpvar(nstre+1,:))
                        valun(ielem)=valun(ielem)/ngaus
                        nullify(lnods,ldofs)
                    endif
                end do

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (nelem1>0)then
                    DO ielgroup = 1,group1(igroup)%nelgroup
                        ielem = group1(igroup)%list(ielgroup)
                        if (jce1(ielem)==0)then
                            lnods=>element1(ielem)%field(1)%lnods_f
                            ldofs=>element1(ielem)%field(1)%ldofs_f
                            eldis = result_zero(ldofs)
                            valun1(ielem)=0.
                            do igaus=1,ngaus
                                cartd=element1(ielem)%egaus(order_int)%cartd(:,:,igaus)
                                do idimn=1,ndimn
                                    stran(idimn)=0.
                                    do inode=1,nnode
                                        stran(idimn)=stran(idimn)+cartd(idimn,inode)*eldis(ndimn*(inode-1)+idimn)
                                    end do
                                end do
                                tstran=sqrt(sum(stran(1:ndimn)**2))
                                valun1(ielem)=valun1(ielem)+tstran
                            end do
                            valun1(ielem)=valun1(ielem)/ngaus
                            nullify(lnods,ldofs)
                        endif
                    end do
                endif
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                deallocate (stran,cartd,eldis)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            endif           !! for(U,CO)
        endif                            !!for appear>0
    end do !! for igroup
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    stranmax=maxval(valun)
    stranmax1=0.
    if(nelem1>0)stranmax1=maxval(valun1)
    print *,'stranmax=',stranmax,'stranmax1=',stranmax1
    !if(stranmax1>stranmax)stranmax=stranmax1

    allocate(rele(nelem))
    rele=0
    do ielem=1,nelem
        if (ice0(ielem)==0) then
            !         if(valun(ielem)>.98*stranmax)rele(ielem)=1   !for rmesh=1
            if(valun(ielem)>.05*stranmax)rele(ielem)=1  !for rmesh=-1
        endif
    end do
    nelc=sum(rele)
    if (nelc>0)then
        allocate(listnelc(nelc))
        nelc=0
        do ielem=1,nelem
            if (rele(ielem)==1)then
                nelc=nelc+1
                listnelc(nelc)=ielem
            endif
        end do
    endif
    deallocate(rele)

    if(rmesh==1)return

    allocate(rele(nelem1))
    rele=0
    !   return
    do ielem=1,nelem1
        if (jce1(ielem)==0)  then
            if(valun1(ielem)>.5*stranmax)rele(ielem)=1
        endif
    end do
    nelc1=sum(rele)
    if (nelc1>0)then
        allocate(listnelc1(nelc1))
        nelc1=0
        do ielem=1,nelem1
            if (rele(ielem)==1)then
                nelc1=nelc1+1
                listnelc1(nelc1)=ielem
            endif
        end do
    endif
    deallocate(rele)

    END SUBROUTINE find_remesh_stran

    subroutine modf_inpwav         !hxl2006 MIF

    character(80)text
    integer(ink) idofn,idofix,iwavcurve,ifixvar,ifixvar0,ilaymif,inpvar
    real   (irk) dfact1,dfact2,frecoord,bfrecoord
    integer(ink),pointer::lnofixb(:),ldofixb(:)

    do idofix=1,ndofix
        ifixvar=prescrib(idofix)%ifixvar
        ifixvar0=prescrib(idofix)%ifixvar0
        bfrecoord=prescrib(idofix)%bfrecoord
        iwavcurve=earthquake_curve_MIF(ifixvar)
        if(iwavcurve==0)cycle
        lnofixb=>prescrib(idofix)%lnofixb
        ldofixb=>prescrib(idofix)%ldofixb
        inpvar=abs(ifixvar0_inpb)
        if (ifixvar0_inpb==ifixvar0)then !底边界
            do ilaymif=1,nlaymif
                tcurves(iwavcurve)%dtbegin=abs(coord(inpvar,lnofixb(ilaymif))-inpcord)/camif !入射位移波传播至当前层的时间
                call dfact_time_curve(ttime)
                inpru(ldofixb(ilaymif))=tcurves(iwavcurve)%dfact !当前层的入射位移波
                tcurves(iwavcurve)%dtbegin=0.0
            enddo
        else
            do ilaymif=1,nlaymif
                tcurves(iwavcurve)%dtbegin=abs(coord(inpvar,lnofixb(ilaymif))-inpcord)/camif
                call dfact_time_curve(ttime)
                dfact1=tcurves(iwavcurve)%dfact !入射波
                tcurves(iwavcurve)%dtbegin=abs((bfrecoord-inpcord)/camif)+abs(coord(inpvar,lnofixb(ilaymif))-bfrecoord)/camif
                call dfact_time_curve(ttime)
                dfact2=tcurves(iwavcurve)%dfact !入射波传播至自由表面后下行反射波
                inpzi(ldofixb(ilaymif))=dfact1+dfact2 !二者的迭加就是自由场
                tcurves(iwavcurve)%dtbegin=0.0
            enddo
        endif
        nullify(lnofixb,ldofixb)
    enddo

    end subroutine modf_inpwav

    subroutine liquifaction_judge  !20231008

    character(1)field1
    integer(ink)igroup,liquj,ij,index,matno,order_int,ngaus,nalfa,ncycl,  &
        ielgroup,ielem,npeak,id
    real(irk) ceqcy(14),bi,bj,bm,bijm,stmax,ratio,dts,dts_1,neqcy,sigma0,  &
        tshear,tyz,tzx,taij,tbij,alfa1,alfa2,alfac,cycl1,cycl2,tstrength
    real(irk),allocatable::qtime(:),tai(:),tbi(:),speak(:),kliqu(:)
    real(irk),pointer::ta(:,:),tb(:,:),alfai(:),cycli(:)
    DATA ceqcy/3.,2.7,2.4,2.05,1.7,1.4,1.2,1.0,.76,.4,.2,.1,.04,.02/


    rewind(lquunit)

    print *,'in liqu_judge,nliqu=',nliqu
    allocate(shear(ngroup),qtime(nliqu))

    DO igroup =1,ngroup
        liquj=  group(igroup)%liquj
        field1= group(igroup)%fieldid(1:1)
        if(appear(igroup)>0.and.field1=='U'.and.liquj==1) then
            allocate(shear(igroup)%stres(nliqu,group(igroup)%nelgroup))
        endif
    end do

    do i0=1,nliqu
        read(lquunit)qtime(i0)
        !  write(7,*)'qtime=',qtime(i0)
        DO igroup =1,ngroup
            liquj=  group(igroup)%liquj
            field1= group(igroup)%fieldid(1:1)
            if(appear(igroup)>0.and.field1=='U'.and.liquj==1) then
                read(lquunit)shear(igroup)%stres(i0,:)
            endif
        end do
    end do
    !!!!!!分组计算周数

    DO igroup =1,ngroup
        liquj=  group(igroup)%liquj
        field1= group(igroup)%fieldid(1:1)
        if(appear(igroup)>0.and.field1=='U'.and.liquj==1) then
            index=group(igroup)%index
            order_int=elkn(index)%el_field(1)%order_intrules(1)
            ngaus=elkn(index)%ggaus(order_int)%ngaus
            matno=group(igroup)%matno
            if(props(matno)%mechanical%solid%jliqu==0)then
                print *,'stop err in material number,  &
                    should include the parameters for liqufaction'
                stop
            endif

            nalfa=props(matno)%mechanical%solid%scycl%nalfa
            ncycl=props(matno)%mechanical%solid%scycl%ncycl
            ta=>props(matno)%mechanical%solid%scycl%ta
            tb=>props(matno)%mechanical%solid%scycl%tb
            alfai=>props(matno)%mechanical%solid%scycl%alfai
            cycli=>props(matno)%mechanical%solid%scycl%cycli
            allocate(tai(ncycl),tbi(ncycl))
            allocate(speak(nliqu),kliqu(group(igroup)%nelgroup))
            DO ielgroup = 1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)

                npeak=0
                !!计算峰值点个数
                do i0=2,nliqu-1
                    bi=shear(igroup)%stres(i0-1,ielgroup)
                    bj=shear(igroup)%stres(i0,ielgroup)
                    bm=shear(igroup)%stres(i0+1,ielgroup)
                    bijm=(bj-bi)*(bj-bm)
                    !			write(7,*)'ig=',igroup,'ie=',ielgroup,'bi,bj,bm=',bi,bj,bm,'bijm=',bijm
                    if(bijm>0)then
                        if((bj>0..and.(bj-bi)>0.).or.(bj<0..and.(bj-bi)<0.))then
                            npeak=npeak+1
                            !				  write(7,*)'npeak=',npeak
                            speak(npeak)=bj
                        endif
                    endif
                end do
                !        write(7,*)'ig=',igroup,'ie=',ielgroup,'npeak=',npeak,'speak=',speak(1:npeak)
                !!end计算峰值点个数
                !!计算等效周数
                neqcy=0
                if(npeak>0)then
                    !write(7,*)'*******************************npeak=*********************',npeak
                    stmax=maxval(abs(shear(igroup)%stres(:,ielgroup)))
                    speak=abs(speak)/stmax

                    do ij=1,npeak
                        do id=1,14
                            dts=1.-(id-1)*.05
                            dts_1=dts-.05
                            if(id==14)dts_1=0.
                            if(speak(ij)<=dts.and.speak(ij)>dts_1)neqcy=neqcy+.5*ceqcy(id)
                        end do
                    end do
                end if
                !	print *,'neqcy=',neqcy
                element(ielem)%neqcy=neqcy
                !!end计算等效周数
                !! take out the initial vertical normal stress and the horizontal shear stress
                !sigma0=sum(element(ielem)%field(1)%STRES0(ndimn,:))/ngaus
                sigma0=sum(element(ielem)%STRES0(ndimn,:))/ngaus
                if(ndimn==2)then
                    !tshear=sum(element(ielem)%field(1)%STRES0(3,:))/ngaus
                    tshear=sum(element(ielem)%STRES0(3,:))/ngaus
                elseif(ndimn==3)then
                    !tyz=sum(element(ielem)%field(1)%STRES0(5,:))/ngaus
                    !tzx=sum(element(ielem)%field(1)%STRES0(6,:))/ngaus
                    tyz=sum(element(ielem)%STRES0(5,:))/ngaus
                    tzx=sum(element(ielem)%STRES0(6,:))/ngaus

                    tshear=sqrt(tyz**2+tzx**2)
                endif
                ratio=abs(tshear/sigma0)
                !	print *,'tshear=',tshear,'ratio=',ratio

                tai=0;tbi=0.
                if(ratio>=alfai(nalfa)) then
                    tai=ta(nalfa,:)
                    tbi=ta(nalfa,:)
                else
                    do i0=1,nalfa-1
                        alfa1=alfai(i0)
                        alfa2=alfai(i0+1)
                        if(ratio>=alfa1.and.ratio<alfa2) then
                            alfac=(ratio-alfa1)/(alfa2-alfa1)
                            tai(:)=ta(i0,:)+alfac*(ta(i0+1,:)-ta(i0,:))
                            tbi(:)=tb(i0,:)+alfac*(tb(i0+1,:)-tb(i0,:))
                            goto 10
                        endif
                    end do
                endif

10              taij=0.;tbij=0.

                if(neqcy<=cycli(1)) then
                    taij=tai(ncycl)
                    tbij=tbi(ncycl)
                elseif(neqcy>=cycli(ncycl)) then
                    taij=tai(ncycl)
                    tbij=tbi(ncycl)
                else
                    do i0=1,ncycl-1
                        cycl1=cycli(i0)
                        cycl2=cycli(i0+1)
                        if(neqcy>=cycl1.and.neqcy<=cycl2) then
                            alfac=(neqcy-cycl1)/(cycl2-cycl1)
                            taij=tai(i0)+alfac*(tai(i0+1)-tai(i0))
                            tbij=tbi(i0)+alfac*(tbi(i0+1)-tbi(i0))
                            goto 20
                        endif
                    end do
                endif

20              continue

                tstrength=tbij+taij*abs(sigma0)
                kliqu(ielgroup)=tstrength/(.65*stmax)
                !  write(7,*)'ielem=',ielem,'stmax=',stmax,'tstrength=',tstrength,'kliqu=',kliqu(ielgroup)
            end do
            DO ielgroup = 1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                write(disunit,30)igroup,ielem,kliqu(ielgroup)
            end do
            deallocate(tai,tbi)
            deallocate(speak,kliqu)
            nullify(ta,tb)
        endif
    end do
    !!!!!!!
    DO igroup =1,ngroup
        liquj=  group(igroup)%liquj
        if(appear(igroup)>0.and.field1=='U'.and.liquj==1) then
            deallocate(shear(igroup)%stres)
        endif
    end do
    deallocate(shear,qtime)

30  format(2i10,e20.5)

    end subroutine liquifaction_judge

    subroutine permdeform_judge !20231008

    character(100) name,material,SPtype,field1
    integer(ink) igroup,index,matno,nstre,nnode,ielgroup,ielem,igaus,ngaus,order_int,nevab,liquj
    real   (irk) pei,totneq,neqcy,theta,steff,smean,varj3,varj2,toct,p0,p,p3,s,q,qf,ROOT3,strain_s,devr,drr, &
        c1,c2,c3,c4,c5,vj3,phi,COHES,vj2,sint3
    real   (irk),allocatable::unitx(:),SGTOT(:),devia(:),rot(:),ps(:),eldis(:),bmatx(:,:),stran(:),stemp(:), &
        stmin(:),sigx(:),strand(:)
    integer(ink),pointer::ldofs(:)
    real   (irk),pointer::cartd(:,:),shape(:),gpcod(:)


    pei= 3.14159 ;       ROOT3=1.73205080757


    allocate(unitx(3*(ndimn-1)),ps(ndimn))
    ps=0. ; totneq=0.

    unitx=1.
    unitx(ndimn+1:3*(ndimn-1))=0.
    write(pmtunit,*)'ttime=',ttime,' istep=',istep
    DO igroup =1,ngroup
        liquj=  group(igroup)%liquj
        index   =group(igroup)%index
        order_int=elkn(index)%el_field(1)%order_intrules(1)
        ngaus=elkn(index)%ggaus(order_int)%ngaus
        matno   =group(igroup)%matno
        name    =props(matno)%name
        field1= group(igroup)%fieldid(1:1)    ! nzw 2016-02-02   field1(1:1) for UW
        material=props(matno)%mechanical%solid%material
        if(appear(igroup)/=1.or.field1/='U'.or.material/='DUNCANCHANG'.or.liquj/=1)cycle !CR?    !psy.2018.09.20
        nstre   =group(igroup)%nstre
        SPtype  =group(igroup)%SPtype
        nnode   =elkn(index)%el_field(1)%nnode_f
        phi  =props(matno)%mechanical%solid%DuncanChang%phi
        COHES=props(matno)%mechanical%solid%DuncanChang%COHES
        p0   =props(matno)%mechanical%solid%DuncanChang%p0
        DO ielgroup=1,group(igroup)%nelgroup
            ielem=group(igroup)%list(ielgroup)
            neqcy=element(ielem)%neqcy
            !		  totneq=totneq+neqcy
            write(pmtunit,*)'igroup,   ielem,     neqcy'
            write(pmtunit,*) igroup,ielem,neqcy
            ldofs=>element(ielem)%field(1)%ldofs_f
            nevab=size(ldofs)      !   nevab=size(element(ielem)%ldofs)    nzw 2016-02-02   for UW
            allocate(eldis(nevab))
            eldis =result_zero(ldofs)
            do igaus=1,ngaus
                allocate(SGTOT(nstre),devia(nstre),rot(ndimn))
                devia=0. ; rot=0.

                !           SGTOT=element(ielem)%field(1)%gpvar(1:nstre,igaus)
                SGTOT=element(ielem)%stres0(1:nstre,igaus)

                CALL INVART(matno,nstre,DEVIA,SGTOT,THETA,STEFF,SMEAN,vj2,vj3,sint3,rot)
                if (ndimn==3)then
                    varj2=((SGTOT(1)-SGTOT(2))**2+(SGTOT(2)-SGTOT(3))**2+(SGTOT(1)-SGTOT(3))**2)/6.  &
                        +SGTOT(4)**2+SGTOT(5)**2+SGTOT(6)**2
                elseif(ndimn==2)then
                    if(nstre==3)then
                        varj2=((SGTOT(1)-SGTOT(2))**2+SGTOT(2)**2+SGTOT(1)**2)/6.+SGTOT(3)**2
                    elseif(nstre==4)then
                        varj2=((SGTOT(1)-SGTOT(2))**2+(SGTOT(2)-SGTOT(4))**2+(SGTOT(1)-SGTOT(4))**2)/6.  &
                            +SGTOT(3)**2
                    endif
                    !                 if (nstre==4)then
                    !                    varj2=((SGTOT(1)-SGTOT(2))**2+SGTOT(2)**2+SGTOT(1)**2)/6.+SGTOT(3)**2
                else
                    stop 'nstre ! permdeform_judge'
                endif
                !            endif
                toct=sqrt(2*varj2/3)
                ps(3)=-(2.*steff/root3*sin(theta+2*pei/3.)+smean)
                ps(2)=-(2.*steff/root3*sin(theta         )+smean)
                ps(1)=-(2.*steff/root3*sin(theta+4*pei/3.)+smean)
                p=ps(3)
                if(p<p0)p=p0
                p3=p
                Q=ps(1)-ps(3)
                Qf=(2.*p3*sinD(phi)+2.*COHES*COSD(phi))/(1-sind(phi))
                S=Q/QF
                if(s>1.)s=1.

                cartd=>element(ielem)%egaus(order_int)%cartd(:,:,igaus)
                shape=>elkn(index)%ggaus(order_int)%shape(:,igaus)
                gpcod=>element(ielem)%egaus(order_int)%gpcod(:,igaus)
                allocate(bmatx(nstre,nevab),stran(nstre),stemp(nstre),stmin(ndimn))
                bmatx=0. ; stran=0. ; stmin=0.
                call gbmat(SPtype,nnode,bmatx,cartd,gpcod,shape)
                stran=matmul(bmatx,eldis)
                stemp=stran
                stemp(ndimn+1:3*(ndimn-1))=.5*stran(ndimn+1:3*(ndimn-1))
                call main_stran_r( stemp, stmin)
                ! if(ndimn==2)strain_s=abs((stmin(1)-stmin(2)))*0.5  !最大剪应变(2D) !zhao
                ! if(ndimn==3)strain_s=abs((stmin(1)-stmin(3)))*0.5  !最大剪应变(3D)
                if(ndimn==2)strain_s=(stmin(1)-stmin(2))!*0.5  !最大剪应变(2D) !zhao
                !if(ndimn==3)strain_s=(stmin(1)-stmin(3))!*0.5  !最大剪应变(3D)
                if(ndimn==3)strain_s=sqrt(((stmin(1)-stmin(2))**2+(stmin(2)-stmin(3))**2+(stmin(3)-stmin(1))**2)*2)/3 !最大动剪应变 yuanli
                !			 if (strain_s<=0) stop
                !			 print *,'srain_s0'
                allocate(sigx(nstre))
                sigx=0.
                sigx=SGTOT
                !sigx(1:ndimn)=sigx(1:ndimn)-p
                p=-smean
                !if(p<.01)p=.01
                !if(p<p0)p=p0 !yuanli20230802
                sigx(1:ndimn)=-1*sigx(1:ndimn)-p
                sigx(ndimn+1:3*(ndimn-1))=-1*sigx(ndimn+1:3*(ndimn-1))*2.
                if(istep==nstepJP)write(pmtunit,'(2i8,6e18.6)')ielem,igaus,sigx
                write(pmtunit,'(2i8,6e18.6)')ielem,igaus,s,strain_s,toct
                deallocate(SGTOT,devia,rot,bmatx,stran,stemp,stmin,sigx)
                nullify(cartd,shape,gpcod)
            enddo !igaus
            nullify(ldofs)
            deallocate(eldis)
        enddo !ielem
    enddo !igroup

    end subroutine permdeform_judge

    subroutine gamamaxupdate !20231008
    integer(ink) igroup,matno,index,ngaus,nstre,nnode,nevab,ielgroup,ielem,igaus
    integer(ink),pointer::ldofs(:)
    real(irk) gamad,gamad0
    real(irk),allocatable::eldis(:),bmatx(:,:),stran(:),stemp(:),stmin(:)
    character(100) name,material,fieldid

    DO igroup =1,ngroup
        if(appear(igroup)>0) then
            fieldid=group(igroup)%fieldid
            if(fieldid(1:1)=='U')then
                matno = group(igroup)%matno
                index = group(igroup)%index
                ngaus = elkn(index)%ggaus(1)%ngaus
                name  = props(matno)%name
                material=props(matno)%mechanical%solid%material
                if(material/='DUNCANCHANG')cycle
                nstre=  group(igroup)%nstre
                nnode = elkn(index)%el_field(1)%nnode_f
                nevab = nnode*group(igroup)%dof(1)%nfdof
                allocate(eldis(nevab),bmatx(nstre,nevab),stran(nstre),stemp(nstre),stmin(ndimn))
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    ldofs => element(ielem)%field(1)%ldofs_f
                    eldis = result_zero(ldofs)
                    do igaus = 1,ngaus
                        bmatx = element(ielem)%field(1)%bmatx(:,:,igaus)
                        stran = matmul(bmatx,eldis)
                        !未处理平面应力情况
                        stemp=stran
                        stemp(ndimn+1:3*(ndimn-1))=.5*stran(ndimn+1:3*(ndimn-1))
                        call main_stran_r( stemp, stmin)
                        if(ndimn==2)gamad=abs((stmin(1)-stmin(2)))   !*0.5  !最大剪应变(2D) !zhao
                        !if(ndimn==3)strain_s=abs((stmin(1)-stmin(3)))   !*0.5  !最大剪应变(3D)
                        if(ndimn==3)gamad=sqrt(((stmin(1)-stmin(2))**2+(stmin(2)-stmin(3))**2+(stmin(3)-stmin(1))**2)*2)/3 !最大动剪应变 yuanli
                        gamad0 = element(ielem)%field(1)%gamamax0(igaus)
                        if(gamad>=gamad0)then
                            element(ielem)%field(1)%gamamax(igaus) = gamad
                            element(ielem)%field(1)%gamamax0(igaus) = gamad
                        else
                            element(ielem)%field(1)%gamamax(igaus) = gamad0
                            element(ielem)%field(1)%gamamax0(igaus) = gamad0
                        endif
                    enddo
                    nullify(ldofs)
                enddo
                deallocate(eldis,bmatx,stran,stemp,stmin)
            endif
        endif
    enddo

    end subroutine gamamaxupdate

    subroutine readgamamax !20231008
    integer(ink) igroup,matno,index,ngaus,ielgroup,ielem,igaus,i0
    character(100) material,fieldid,text
    Tgamamax0=0
    if(gamamax==1)then
        return
    elseif(gamamax==2)then
        DO igroup =1,ngroup
            if(appear(igroup)>0) then
                fieldid=group(igroup)%fieldid
                if(fieldid(1:1)=='U')then
                    matno = group(igroup)%matno
                    index = group(igroup)%index
                    ngaus = elkn(index)%ggaus(1)%ngaus
                    material=props(matno)%mechanical%solid%material
                    if(material/='DUNCANCHANG')cycle
                    read(gamamaxunit,*)i0
                    DO ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        do igaus = 1,ngaus
                            read(gamamaxunit,*)i0,i0,element(ielem)%field(1)%gamamax_ini(igaus),element(ielem)%field(1)%gamamax_error(igaus)
                        enddo
                    enddo
                endif
            endif
        enddo
        read(gamamaxunit,*) text,i0,text,Tgamamax0,text,i0
    endif


    end subroutine readgamamax

    subroutine writegamamax !20231008
    integer(ink) igroup,matno,index,ngaus,ielgroup,ielem,igaus
    real(irk) gamamax_error_max,gamamax_ratio
    character(100) material,fieldid
    rewind(gamamaxunit)
    Tgamamax=0
    DO igroup =1,ngroup
        if(appear(igroup)>0) then
            fieldid=group(igroup)%fieldid
            if(fieldid(1:1)=='U')then
                matno = group(igroup)%matno
                index = group(igroup)%index
                ngaus = elkn(index)%ggaus(1)%ngaus
                material=props(matno)%mechanical%solid%material
                if(material/='DUNCANCHANG')cycle
                write(gamamaxunit,*)igroup
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    do igaus = 1,ngaus
                        element(ielem)%field(1)%gamamax_error(igaus) = abs((element(ielem)%field(1)%gamamax(igaus)-element(ielem)%field(1)%gamamax_ini(igaus))/element(ielem)%field(1)%gamamax(igaus))
                        element(ielem)%field(1)%gamamax_ini(igaus) = element(ielem)%field(1)%gamamax(igaus)
                        Tgamamax=Tgamamax+element(ielem)%field(1)%gamamax(igaus)
                        write(gamamaxunit,*)ielem,igaus,element(ielem)%field(1)%gamamax_ini(igaus),element(ielem)%field(1)%gamamax_error(igaus)
                        if(gamamax_error_max<element(ielem)%field(1)%gamamax_error(igaus))gamamax_error_max = element(ielem)%field(1)%gamamax_error(igaus)
                    enddo
                enddo
            endif
        endif
    enddo
    gamamax_ratio=abs(Tgamamax-Tgamamax0)/Tgamamax
    write(gamamaxunit,*) '结点最大相对误差',gamamax_error_max,'总误差',Tgamamax,'总相对误差',gamamax_ratio
    close(gamamaxunit)
    end subroutine writegamamax

    SUBROUTINE inivdval !20220721
    character(10) fieldid
    integer(ink) igroup,ielgroup,ielem,matno,order_int,index,igaus,ngaus

    DO igroup =1,ngroup

        if(appear(igroup)>0) then

            fieldid=group(igroup)%fieldid

            if(fieldid(1:2)=='UW'.or.fieldid(1:1)=='U') then		!nzw 2006-06-19 for PZ to U field
                matno = group(igroup)%matno
                index = group(igroup)%index

                if(props(matno)%mechanical%solid%material=='SandPZ'.or. &
                    props(matno)%mechanical%solid%material=='ClayPZ')then

                    order_int=elkn(index)%el_field(1)%order_intrules(1)

                    ngaus = elkn(index)%ggaus(order_int)%ngaus

                    DO ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        element(ielem)%egaus(order_int)%iload0=0
                        element(ielem)%egaus(order_int)%vdval0=0.
                        element(ielem)%egaus(order_int)%iload =0
                        element(ielem)%egaus(order_int)%vdval =0.


                        if(ndimn==2)then

                            do igaus=1,ngaus

                                element(ielem)%egaus(order_int)%vdval0(5,igaus)=           &
                                    -(element(ielem)%field(1)%gpvar0(1,igaus)+         &
                                    element(ielem)%field(1)%gpvar0(2,igaus)+           &
                                    element(ielem)%field(1)%gpvar0(4,igaus))/3.

                                element(ielem)%egaus(order_int)%vdval(5,igaus)=            &
                                    -(element(ielem)%field(1)%gpvar0(1,igaus)+         &
                                    element(ielem)%field(1)%gpvar0(2,igaus)+           &
                                    element(ielem)%field(1)%gpvar0(4,igaus))/3.
                            enddo

                        elseif(ndimn==3)then
                            do igaus=1,ngaus

                                element(ielem)%egaus(order_int)%vdval0(5,igaus)=           &
                                    -(element(ielem)%field(1)%gpvar0(1,igaus)+         &
                                    element(ielem)%field(1)%gpvar0(2,igaus)+           &
                                    element(ielem)%field(1)%gpvar0(3,igaus))/3.

                                element(ielem)%egaus(order_int)%vdval(5,igaus)=            &
                                    -(element(ielem)%field(1)%gpvar0(1,igaus)+         &
                                    element(ielem)%field(1)%gpvar0(2,igaus)+           &
                                    element(ielem)%field(1)%gpvar0(3,igaus))/3.
                            enddo
                        else
                            write(*,*) 'in sub inivaval, ndimn should be 2 or 3'
                            stop
                        endif

                        if(props(matno)%mechanical%solid%material=='ClayPZ')then
                            element(ielem)%egaus(order_int)%vdval0(3,:)=props(matno)%mechanical%solid%ClayPZ%d(10)
                            element(ielem)%egaus(order_int)%vdval(3,:) =props(matno)%mechanical%solid%ClayPZ%d(10)
                            element(ielem)%egaus(order_int)%vdval0(4,:)=props(matno)%mechanical%solid%ClayPZ%d(10)
                            element(ielem)%egaus(order_int)%vdval(4,:) =props(matno)%mechanical%solid%ClayPZ%d(10)
                        elseif(props(matno)%mechanical%solid%material=='SandPZ')then
                            element(ielem)%egaus(order_int)%vdval0(2,:)=props(matno)%mechanical%solid%SandPZ%d(12)   !Hu0
                            element(ielem)%egaus(order_int)%vdval (2,:)=props(matno)%mechanical%solid%SandPZ%d(12)
                        endif
                    end do
                endif
            endif

        endif
    end do

    END  SUBROUTINE inivdval  !20220721

    subroutine strain_for_steel_bar
    integer(ink) ielem,index,nnode,igroup,nevab,ipoin
    real   (irk) dl
    integer(ink),pointer::lnods(:),ldofs(:)
    real   (irk),pointer::rotation(:,:)
    real   (irk),allocatable::eldis(:)

    pstrain=0.

    do ielem=1,nelem
        index =element(ielem)%index
        if (index/=20.and.index/=1)cycle  !20200628
        nnode =elkn(index)%nnode
        igroup=element(ielem)%group
        if(appear(igroup)==0)cycle
        lnods=>element(ielem)%field(1)%lnods_f
        ldofs=> element(ielem)%field(1)%ldofs_f
        rotation=>element(ielem)%rotation
        nevab=size(ldofs) ; allocate(eldis(nevab)) ; eldis=0.
        dl=sqrt(sum((coord(:,lnods(2))-coord(:,lnods(1)))**2))
        eldis=result_zero(ldofs)
        if(nlocalbeam==0)then  !20230906
            !pstrain(lnods)=pstrain(lnods)+((eldis(nevab/2+1:nevab/2+ndimn)-eldis(1:ndimn)).d.rotation(1,:))/dl
            pstrain(lnods)=pstrain(lnods)+((eldis(ndimn+1:ndimn)-eldis(1:ndimn)).d.rotation(1,:))/dl !20230906
        else
            pstrain(lnods)=pstrain(lnods)+(eldis(nevab/2+1)-eldis(1))/dl
        endif
        nullify(lnods,ldofs,rotation)
        deallocate(eldis)
    enddo
    do ipoin=1,npoin
        if(icpspring(ipoin)<1)cycle
        pstrain(ipoin)=pstrain(ipoin)/real(icpspring(ipoin))
    enddo

    end subroutine strain_for_steel_bar

    subroutine change_list(listnode,index,nnode)
    integer(ink) index,nnode,listnode(nnode),nodeface(4),listnode0(nnode)
    integer(ink),allocatable::ipi(:)
    allocate(ipi(npoin))
    ipi=0
    ipi(listnode)=1
    select case(index)
    case(5,22)
        if(sum(ipi)==3)call change4(listnode)   !对退化的四边形处理
    case(9)
        if(sum(ipi)==6)then  !对退化的六面体处理
            ipi=0
            ipi(listnode(1:4))=1

            if(sum(ipi)==3)then
                nodeface=listnode(1:4)
                call change4(nodeface)
                listnode(1:4)=nodeface
                nodeface=listnode(5:8)
                call change4(nodeface)
                listnode(5:8)=nodeface
                deallocate(ipi)
                return
            endif

            ipi=0
            ipi(listnode(1))=1;ipi(listnode(5))=1;ipi(listnode(6))=1;ipi(listnode(2))=1
            if(sum(ipi)==3)then
                listnode0=listnode
                listnode(1)=listnode0(1)
                listnode(2)=listnode0(5)
                listnode(3)=listnode0(6)
                listnode(4)=listnode0(2)
                listnode(5)=listnode0(4)
                listnode(6)=listnode0(8)
                listnode(7)=listnode0(7)
                listnode(8)=listnode0(3)

                nodeface=listnode(1:4)
                call change4(nodeface)
                listnode(1:4)=nodeface
                nodeface=listnode(5:8)
                call change4(nodeface)
                listnode(5:8)=nodeface
                deallocate(ipi)
                return
            endif

            ipi=0
            ipi(listnode(2))=1;ipi(listnode(6))=1;ipi(listnode(7))=1;ipi(listnode(3))=1
            if(sum(ipi)==3)then
                listnode0=listnode
                listnode(1)=listnode0(2)
                listnode(2)=listnode0(6)
                listnode(3)=listnode0(7)
                listnode(4)=listnode0(3)
                listnode(5)=listnode0(1)
                listnode(6)=listnode0(5)
                listnode(7)=listnode0(8)
                listnode(8)=listnode0(4)

                nodeface=listnode(1:4)
                call change4(nodeface)
                listnode(1:4)=nodeface
                nodeface=listnode(5:8)
                call change4(nodeface)
                listnode(5:8)=nodeface
                deallocate(ipi)
                return
            endif
        endif
    end select
    if(allocated(ipi))deallocate(ipi)
    end subroutine change_list

    subroutine change4(listnode)
    integer(ink) listnode(4),ipoin

    if(listnode(1)==listnode(2))then
        ipoin=listnode(1)
        listnode(1:2)=listnode(3:4)
        listnode(3:4)=ipoin
    elseif(listnode(2)==listnode(3))then
        ipoin=listnode(1)
        listnode(1)=listnode(4)
        listnode(2)=ipoin
        listnode(4)=listnode(3)
    elseif(listnode(1)==listnode(4))then
        ipoin=listnode(1)
        listnode(1:2)=listnode(2:3)
        listnode(3)=listnode(4)
    endif

    end subroutine change4

    subroutine GHM2ADINA
    integer(ink),allocatable::ipiface(:),list_fix(:),val_fix(:),fix_xyz(:)
    integer(ink),pointer::lnods(:),nlist(:)
    real(irk),pointer::pxyz(:),fxyz(:,:)
    character*20 gname,name,material,criteria,text,sptype,type_curve
    integer(ink) adnunit,imat,nnode,nset,iplgroup,iedge,inode,edimn,ndofn,iset,   &
        ifixsets,nfixsets,ifixnods,nfixnods,pset,ntime,itime,nload_mass,iload_mass,idimn
    real(irk) density,e,nu,alpha,frict_angle,sigma0,pvalue,rot(3),dp_alfa,dp_k,r0,r1,deltatime

    adnunit=101
    open(adnunit,  file=trim(probn)//'.in')

    write(adnunit,1000)'DATABASE NEW SAVE=NO PROMPT=NO'
    write(adnunit,1000)'FEPROGRAM ADINA'
    write(adnunit,1000)'CONTROL FILEVERSION=V85'

    if(type_problem=='Q')write(adnunit,'(a)')'MASTER ANALYSIS=STATIC MODEX=EXECUTE TSTART=0.00000000000000,'
    if(type_problem=='F')write(adnunit,'(a)')'MASTER ANALYSIS=DYNAMIC-DIRECT-INTEGRATION MODEX=EXECUTE TSTART=0.00000000000000,'

    if(ndimn==2.and.MDOFN==2)then
        write(adnunit,'(a)')'     IDOF=100111 OVALIZAT=NONE FLUIDPOT=AUTOMATIC CYCLICPA=1,'
    elseif(ndimn==3.and.MDOFN==3)then
        write(adnunit,'(a)')'     IDOF=111 OVALIZAT=NONE FLUIDPOT=AUTOMATIC CYCLICPA=1,'
    endif
    write(adnunit,'(a)')'     IPOSIT=STOP REACTION=YES INITIALS=NO FSINTERA=NO IRINT=DEFAULT,'
    write(adnunit,'(a)')'     CMASS=NO SHELLNDO=AUTOMATIC AUTOMATI=OFF SOLVER=SPARSE,'
    write(adnunit,'(a)')'     CONTACT-=CONSTRAINT-FUNCTION TRELEASE=0.00000000000000,'
    write(adnunit,'(a)')'     RESTART-=NO FRACTURE=NO LOAD-CAS=NO LOAD-PEN=NO MAXSOLME=0,'
    write(adnunit,'(a)')'     MTOTM=2 RECL=3000 SINGULAR=YES STIFFNES=0.000100000000000000,'
    write(adnunit,'(a)')"     MAP-OUTP=NONE MAP-FORM=NO NODAL-DE='' POROUS-C=NO ADAPTIVE=0,"
    write(adnunit,'(a)')'     ZOOM-LAB=1 AXIS-CYC=0 PERIODIC=NO VECTOR-S=GEOMETRY EPSI-FIR=NO,'
    write(adnunit,'(a)')'     STABILIZ=NO STABFACT=1.00000000000000E-12 RESULTS=PORTHOLE,'
    write(adnunit,'(a)')'     FEFCORR=NO BOLTSTEP=1 EXTEND-S=YES CONVERT-=NO DEGEN=YES'

    !write material information
    write(adnunit,1000)'* define material proterties'
    do imat=1,nmats
        if(associated(props(imat)%mechanical%solid))then
            material=props(imat)%mechanical%solid%material
            density=props(imat)%mechanical%solid%density
            if(Bparameter/=0.and.props(imat)%mechanical%solid%ie/=0)then !20190810
                e=xvalue(props(imat)%mechanical%solid%ie)
            else
                e=props(imat)%mechanical%solid%e !exx !
            endif
            if(Bparameter/=0.and.props(imat)%mechanical%solid%iNu/=0)then
                Nu=xvalue(props(imat)%mechanical%solid%iNu)
            else
                Nu=props(imat)%mechanical%solid%Nu !uxx !
            endif
            alpha=props(imat)%mechanical%solid%alfa
            select    case(material)
            case('ELASTIC_ISOTROPIC')
                write(adnunit,'(a,i3,a,e16.6,a,f5.3,a)')'MATERIAL ELASTIC NAME=',imat,' E=',e,' NU=',nu,','
                write(adnunit,'(2(a,e16.6),a)')"     DENSITY=",density," ALPHA=",alpha," MDESCRIP='NONE'"
            case('CLASSICALEP')
                criteria=props(imat)%mechanical%solid%ClassicalEP%criteria
                if(criteria(1:2)=='MC')then
                    frict_angle=props(imat)%mechanical%solid%ClassicalEP%frict_angle
                    sigma0=props(imat)%mechanical%solid%ClassicalEP%sigma0
                    write(adnunit,'(a,i3,a,e16.6,a)')'MATERIAL MOHR-COULOMB NAME=',imat,' E=',e,','
                    write(adnunit,'(a,f5.3,a,f8.3,a)')'     NU=',nu,' PHI=',frict_angle,' PSI=0.00000000000000,'
                    write(adnunit,'(a,e16.6,a)')'     COH=',sigma0,' TCUT=0.00000000,'
                    write(adnunit,'(a,e16.6,a)')"     DENSITY=",density," DILATION=NO MDESCRIP='NONE'"
                elseif(criteria(1:2)=='DP')then
                    frict_angle=props(imat)%mechanical%solid%ClassicalEP%frict_angle
                    sigma0=props(imat)%mechanical%solid%ClassicalEP%sigma0
                    !if(criteria=='DP1')then !外顶点
                    !    dp_alfa=2.0*sind(frict_angle)/(sqrt(3.0)*(3.0-sind(frict_angle)))
                    !    dp_k=   6.0*sigma0*cosd(frict_angle)/(sqrt(3.0)*(3.0-sind(frict_angle)))
                    !elseif(criteria=='DP2')then !内顶点  DP3内切
                    dp_alfa=2.0*sind(frict_angle)/(sqrt(3.0)*(3.0+sind(frict_angle)))
                    dp_k=   6.0*sigma0*cosd(frict_angle)/(sqrt(3.0)*(3.0+sind(frict_angle)))
                    !else
                    !    write(*,*)'**********ERROR************'
                    !    write(*,*)'GHM2ADINA. Criteria should be DP1 or DP2 for matno=',imat
                    !    write(*,*)'**********ERROR************'
                    !    stop
                    !endif
                    write(adnunit,'(a,i3,a,e16.6,a)')'MATERIAL DRUCKER-PRAGER NAME=',imat,' E=',e,','
                    write(adnunit,'(a,f5.3,a,e16.6,a)')'     NU=',nu,' ALPHA=',dp_alfa,','
                    write(adnunit,'(a,e16.6,a)')'     KYIELD=',dp_k,' WCAP=-0.100000000000000,'
                    write(adnunit,'(a)')'     DCAP=-0.100000000000000 TCUT=100000.000000000,'
                    write(adnunit,'(a)')'     ICPOS=0.00000000000000 RCAP=0.00000000000000,'
                    write(adnunit,'(a,e16.6,a)')'     DENSITY=',density,' BETA=0.00000000000000 POTENTIA=NO,'
                    write(adnunit,'(a)')"     MDESCRIP='NONE'"
                endif
                case default
                write(adnunit,1000)'*This material is not included in adina. matno=',imat
                write(*,1000)'*This material is not included in adina. matno=',imat
                stop
            end select

        elseif(associated(props(imat)%mechanical%fluid))then
        else
        endif

        write(adnunit,1000)
    end do
    !write coordinate information
    write(adnunit,1000)'COORDINATES NODE SYSTEM=0'
    write(adnunit,1000)'@CLEAR'
    do ipoin=1,npoin
        if(ndimn==2)write(adnunit,1001)ipoin,0.0,coord(1:ndimn,ipoin)    !convert to YZ plane
        if(ndimn==3)write(adnunit,1001)ipoin,coord(1:ndimn,ipoin)
    end do
    write(adnunit,'(a)')'@'
    write(adnunit,'(a)')'*'

    !write element information
    write(adnunit,1000)'* define element information'

    do igroup=1,ngroup
        name=group(igroup)%kname
        index=group(igroup)%index
        nnode=elkn(index)%nnode
        sptype=group(igroup)%sptype
        if(sptype=='PS')sptype='STRESS'
        if(sptype=='PE')sptype='STRAIN'
        matno=matno_process(igroup,iblks)
        nnode=elkn(index)%nnode
        if(index==1.or.index==2)then
            write(adnunit,1002)
        elseif(index==19.or.index==20.or.index==21)then
            write(adnunit,1002)
        elseif(index==3.or.index==4)then
            write(adnunit,1002)
        elseif(index==5.or.index==6.or.index==12.or.index==16)then
            write(text,'(i4)')igroup
            write(adnunit,'(8a)')'EGROUP TWODSOLID NAME=',trim(ADJUSTL(text)),' SUBTYPE=',trim(sptype),' DISPLACE=DEFAULT,'
            write(adnunit,'(a,i3,a)')'     STRAINS=DEFAULT MATERIAL=',matno,' INT=DEFAULT RESULTS=STRESSES,'
            write(adnunit,'(a)')'     DEGEN=YES FORMULAT=0 STRESSRE=GLOBAL INITIALS=NONE FRACTUR=NO,'
            write(adnunit,'(a)')'     CMASS=DEFAULT STRAIN-F=0 UL-FORMU=DEFAULT PNTGPS=0 NODGPS=0,'
            write(adnunit,'(a)')'     LVUS1=0 LVUS2=0 SED=NO RUPTURE=ADINA INCOMPAT=DEFAULT,'
            write(adnunit,'(a)')'     TIME-OFF=0.00000000000000 POROUS=NO WTMC=1.00000000000000,'
            write(adnunit,'(a)')"     OPTION=NONE DESCRIPT='NONE' THICKNES=1.00000000000000,"
            write(adnunit,'(a)')'     PRINT=DEFAULT SAVE=DEFAULT TBIRTH=0.00000000000000,'
            write(adnunit,'(a)')'     TDEATH=0.00000000000000'
        elseif(index==22.or.index==26)then  !20230910
            write(adnunit,1002)
        elseif(index==7.or.index==8.or.index==13.or.index==17)then
            write(adnunit,1002)
        elseif(index==9.or.index==10.or.index==14.or.index==18)then
            write(text,'(i4)')igroup
            write(adnunit,'(3a,i3,a)')'EGROUP THREEDSOLID NAME=',trim(ADJUSTL(text)),' DISPLACE=DEFAULT STRAINS=DEFAULT MATERIAL=',matno,','
            write(adnunit,'(a)')'     RSINT=DEFAULT TINT=DEFAULT RESULTS=STRESSES DEGEN=YES FORMULAT=0,'
            write(adnunit,'(a)')'     STRESSRE=GLOBAL INITIALS=NONE FRACTUR=NO CMASS=DEFAULT,'
            write(adnunit,'(a)')'     STRAIN-F=0 UL-FORMU=DEFAULT LVUS1=0 LVUS2=0 SED=NO RUPTURE=ADINA,'
            write(adnunit,'(a)')'     INCOMPAT=DEFAULT TIME-OFF=0.00000000000000 POROUS=NO,'
            write(adnunit,'(a)')"     WTMC=1.00000000000000 OPTION=NONE DESCRIPT='NONE' PRINT=DEFAULT,"
            write(adnunit,'(a)')'     SAVE=DEFAULT TBIRTH=0.00000000000000 TDEATH=0.00000000000000'
        endif

        write(adnunit,'(3a)')'ENODES SUBSTRUC=0 GROUP=',trim(ADJUSTL(text)),' NNODES=32'
        write(adnunit,'(a)')'@CLEAR'

        DO ielgroup = 1,group(igroup)%nelgroup
            ielem = group(igroup)%list(ielgroup)
            lnods=>element(ielem)%field(1)%lnods_f
            call change_list(lnods,index,nnode)
            if(index==1)then
            elseif(index==5)then
                write(adnunit,'(5i8,a)')ielem,lnods,' 0 0 0 0 0'
            elseif(index==9)then
                write(adnunit,'(9i8,a)')ielem,lnods,' 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
            endif
            nullify(lnods)
        end do
    end do

    !write time curve
    do itcurve=1,ntcurve
        type_curve=tcurves(itcurve)%type_curve
        ntime=tcurves(itcurve)%ntime
        if(type_curve=='LINEAR'.or.type_curve=='SEISMIC')then
            write(adnunit,'(a,i2,a)')'TIMEFUNCTION NAME=',itcurve,' IFLIB=1 FPAR1=0.0,FPAR2=0.0 FPAR3=0.0,FPAR4=0.0 FPAR5=0.0,FPAR6=0.0'
            write(adnunit,'(a)')'@CLEAR'
            do itime=1,ntime
                if(type_curve=='SEISMIC')write(adnunit,'(2e16.6)')itime*tcurves(itcurve)%dtrec,tcurves(itcurve)%dfact_curve(itime)*tcurves(itcurve)%ample
                if(type_curve=='LINEAR')write(adnunit,'(2e16.6)')tcurves(itcurve)%ttime_curve(itime),tcurves(itcurve)%dfact_curve(itime)
            end do
            write(adnunit,'(a)')'@'
        endif
    end do

    !write node load
    nset=0
    do iplgroup=1,nplgroup

        pxyz=>pload(iplgroup)%pxyz
        rot=0.
        pvalue=sqrt(dot_product(pxyz,pxyz))
        if(pvalue>=1e-6)then
            nset=nset+1
            !define nole set
            nlist=>pload(iplgroup)%list
            write(adnunit,'(a,i3,a)')"NODESET NAME=",nset," ALL-EXT=NO DESCRIPT='NONE' OPTION=NODE GROUP=0 ZONE='',ELSET=0 TARGET=0"
            write(adnunit,'(a)')'@CLEAR'
            do inode=1,size(nlist,1)
                write(adnunit,'(i8,2i3)')nlist(inode),0,1
            end do
            write(adnunit,'(a)')'@'
            nullify(nlist)

            !define node load value and vector
            rot(1:ndimn)=pxyz/pvalue
            if(ndimn==2)then
                rot(3)=rot(2)
                rot(2)=rot(1)
                rot(1)=0.
            endif
            write(adnunit,'(a,i3,a,e16.6,a,f8.5,a)')'LOAD FORCE NAME=',iplgroup,' MAGNITUD=',pvalue,' FX=',rot(1),','
            write(adnunit,'(2(a,f8.5))')'     FY=',rot(2),' FZ=',rot(3)
        endif
        nullify(pxyz)
    end do

    !wirte face load
    allocate(ipiface(npoin),fxyz(ndimn,npoin))
    ipiface=0;fxyz=0.
    do iedge=1,nedge
        nnode=edges(iedge)%nnode
        index=edges(iedge)%index
        edimn=elkn(index)%ndimn
        nlist=>edges(iedge)%lnode
        ipiface(nlist)=1
        ndofn=ndimn
        do inode=1,nnode
            idofn=(inode-1)*ndofn
            fxyz(1:ndimn,nlist(inode))=fxyz(1:ndimn,nlist(inode))+edgeload(iedge)%edload(idofn+1:idofn+edimn+1)
        end do
        nullify(nlist)
    end do

    do ipoin=1,npoin
        if(ipiface(ipoin)==0)cycle
        rot=0.
        pvalue=sqrt(dot_product(fxyz(1:ndimn,ipoin),fxyz(1:ndimn,ipoin)))
        if(pvalue>=1e-6)then
            !define nole set
            nset=nset+1
            write(adnunit,'(a,i3,a)')"NODESET NAME=",nset," ALL-EXT=NO DESCRIPT='NONE' OPTION=NODE GROUP=0 ZONE='',"
            write(adnunit,'(a)')'ELSET=0 TARGET=0'
            write(adnunit,'(a)')'@CLEAR'
            write(adnunit,'(i8,2i3)')ipoin,0,1
            write(adnunit,'(a)')'@'

            rot(1:ndimn)=fxyz(1:ndimn,ipoin)/pvalue
            if(ndimn==2)then
                rot(3)=rot(2)
                rot(2)=rot(1)
                rot(1)=0.
            endif

            !define node load value and vector
            write(adnunit,'(a,i3,a,e16.6,a,f8.5,a)')'LOAD FORCE NAME=',nset,' MAGNITUD=',pvalue,' FX=',rot(1),','
            write(adnunit,'(2(a,f8.5))')'     FY=',rot(2),' FZ=',rot(3)
        endif
    end do

    if(.not.allocated(earthquake_curve))then
        allocate(earthquake_curve(ndimn))
        earthquake_curve=0
    endif
    read(mainunit,*)text
    !read(mainunit,*)nincs,r0,r1,earthquake_curve
    read(mainunit,*)nincs,earthquake_curve
    !define body load
    nload_mass=0
    if(any(tcurvegravity/=0)) then
        nload_mass=nload_mass+1
        write(adnunit,'(a,i2,a)')'LOAD MASS-PROPORTIONAL NAME=',nload_mass,' MAGNITUD=9.81,AX=0.0 AY=0.0 AZ=-1.0,INTERPRE=BODY-FORCE'
    endif
    if(earthquake_curve(1)/=0)then
        nload_mass=nload_mass+1
        write(adnunit,'(a,i2,a)')'LOAD MASS-PROPORTIONAL NAME=',nload_mass,' MAGNITUD=1.0,AX=1.0 AY=0.0 AZ=0.0,INTERPRE=GROUND-ACCELERATION'
    endif
    if(earthquake_curve(2)/=0)then
        nload_mass=nload_mass+1
        write(adnunit,'(a,i2,a)')'LOAD MASS-PROPORTIONAL NAME=',nload_mass,' MAGNITUD=1.0,AX=0.0 AY=1.0 AZ=0.0,INTERPRE=GROUND-ACCELERATION'
    endif
    if(ndimn==3)then
        if(earthquake_curve(3)/=0)then
            nload_mass=nload_mass+1
            write(adnunit,'(a,i2,a)')'LOAD MASS-PROPORTIONAL NAME=',nload_mass,' MAGNITUD=1.0,AX=0.0 AY=0.0 AZ=1.0,INTERPRE=GROUND-ACCELERATION'
        endif
    endif

    !apply node load to node set
    write(adnunit,'(a)')'APPLY-LOAD BODY=0'
    write(adnunit,'(a)')'@CLEAR'
    do iset=1,nset
        write(adnunit,'(i3,a,i3,a,i3,a)')iset,"  'FORCE' ",iset,"  'NODE' ",iset," 0 1 0.0 0 -1 0 0 0  'NO',0.0 0.0 1 0"
    end do
    nload_mass=0
    if(any(tcurvegravity/=0)) then
        nload_mass=nload_mass+1
        write(adnunit,'(i3,a,i3,a,i3,a)')iset,"  'MASS-PROPORTIONAL' ",nload_mass,"  'MODEL' 0 0 1 0.0 0 -1 0 0 0,'NO' 0.0 0.0 1 0"
    endif
    do idimn=1,ndimn
        if(earthquake_curve(idimn)==0)cycle
        nload_mass=nload_mass+1
        iset=iset+1
        write(adnunit,'(i3,a,i3,a,i3,a)')iset,"  'MASS-PROPORTIONAL' ",nload_mass,"  'MODEL' 0 0 ",earthquake_curve(idimn)," 0.0 0 -1 0 0 0,'NO' 0.0 0.0 1 0"
    end do
    write(adnunit,'(a)')'@'
    deallocate(ipiface,fxyz,earthquake_curve)

    !wirte constrains
    rewind(punit)
    pset=nset
    read(punit,*)text
    read(punit,*)nfixsets
    allocate(fix_xyz(nfixsets))
    fix_xyz=0
    do ifixsets=1,nfixsets
        read(punit,*)fix_xyz(ifixsets),nfixnods
        if(nfixnods<=0)cycle
        allocate(list_fix(nfixnods),val_fix(nfixnods))
        list_fix=0.;val_fix=0.
        read(punit,*)list_fix(1:nfixnods)
        read(punit,*)val_fix(1:nfixnods)
        pset=pset+1
        !define node set
        write(adnunit,'(a,i3,a)')"NODESET NAME=",pset," ALL-EXT=NO DESCRIPT='NONE' OPTION=NODE GROUP=0 ZONE='',ELSET=0 TARGET=0"
        write(adnunit,'(a)')'@CLEAR'
        do ifixnods=1,nfixnods
            write(adnunit,'(i8,2i3)')list_fix(ifixnods),0,1
        end do
        write(adnunit,'(a)')'@'
        deallocate(list_fix,val_fix)
    end do

    if(ndimn==3)then
        write(adnunit,'(a)')'FIXITY NAME=XD'
        write(adnunit,'(a)')'@CLEAR'
        write(adnunit,'(a)')" 'X-TRANSLATION' 'OVALIZATION'"
        write(adnunit,'(a)')'@'
    endif
    write(adnunit,'(a)')'FIXITY NAME=YD'
    write(adnunit,'(a)')'@CLEAR'
    write(adnunit,'(a)')" 'Y-TRANSLATION' 'OVALIZATION'"
    write(adnunit,'(a)')'@'
    write(adnunit,'(a)')'FIXITY NAME=ZD'
    write(adnunit,'(a)')'@CLEAR'
    write(adnunit,'(a)')" 'Z-TRANSLATION' 'OVALIZATION'"
    write(adnunit,'(a)')'@'


    write(adnunit,'(a)')'*'
    write(adnunit,'(a)')'FIXBOUNDARY NODE-SET FIXITY=ALL'
    write(adnunit,'(a)')'@CLEAR'
    do iset=nset+1,pset
        if(ndimn==2)then
            if(fix_xyz(iset-nset)==1)write(adnunit,'(i8,a)')iset,"  'YD'"
            if(fix_xyz(iset-nset)==2)write(adnunit,'(i8,a)')iset,"  'ZD'"
        else
            if(fix_xyz(iset-nset)==1)write(adnunit,'(i8,a)')iset,"  'XD'"
            if(fix_xyz(iset-nset)==2)write(adnunit,'(i8,a)')iset,"  'YD'"
            if(fix_xyz(iset-nset)==3)write(adnunit,'(i8,a)')iset,"  'ZD'"
        endif
    end do
    write(adnunit,'(a)')'@'

    !write load step
    do iincs=1,nincs
        read(mainunit,*)iset,deltatime,iset,iset,nstep
        read(mainunit,*) !text
        if(type_problem=='F')read(mainunit,*) !text
        write(adnunit,'(a)')'*'
        write(adnunit,'(a)')'TIMESTEP NAME=DEFAULT'
        write(adnunit,'(a)')'@CLEAR'
        write(adnunit,'(i5,e16.6)')nstep,deltatime
        write(adnunit,'(a)')'@'
    end do
1000 format(a,i5,e15.6)
1001 format(i8,3e15.6)
1002 format(a,i5,a)
1003 format(a,10(',',i8))
    close(adnunit)
    end subroutine GHM2ADINA
