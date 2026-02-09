    subroutine process_analysis   !20190810
    implicit none
    integer            igapb,igaps,ipairs,npairs,igapbf,npgblock,i0,ij,ipoin


    ttime=lttime
    print *,'lblks=',lblks,'lincs=',lincs,'ttime=',ttime
    if(Bparameter/=0)then
        result_zero=0.
        if(allocated(result_first))result_first=0.0
        if(allocated(result_second))result_second=0.0
        if(allocated(torel))               torel=0.0
        if(allocated(toforl))              toforl=0.0
        if(nflow/=0)then
            allocate(flowrate(npoin))
            flowrate=0.
        endif


        do igapb=1,ngapb   !tcl 2009/10/11
            if(gapb(igapb)%nrdof==0)cycle
            gapb(igapb)%rdisp_inc=0.
            gapb(igapb)%rdisp_zero=0.
        enddo

        do igaps=1,ngaps
            npairs=gaps(igaps)%npairs
            do ipairs=1,npairs
                if(gaps(igaps)%pair_process(ipairs)==0)cycle  !20200331
                gaps(igaps)%dxyz0(:,ipairs)=0.
                gaps(igaps)%dxyz(:,ipairs)=0.
                gaps(igaps)%ctforce0(:,ipairs)=0.
                gaps(igaps)%ctforce(:,ipairs)=0.
            end do
        end do

        call gpvar_initial  !201605
    endif

    call external_load_1
    do iblks=lblks+1,runblks
        write(7,*)'iblks=',iblks,'runblks=',runblks
        if(allocated(torel))               torel=0.0  !20201121

        if(outplot=='GIDL')call OUT_GID_BIN_MESH

        !if(iblks==uwcpl.or.iblks==nblks)
        ! call gpvar_change  !zhao 2007.04.10  !2016/04/25

        mdiv=1
        print *,'iblks=',iblks
        if(type_load=='DISCONTROL'.and.iblks==1)xload=0.
        if(type_load=='DISCONTROL'.and.iblks>1) xload=yload
        write(chkunit,*)'iblks=',iblks


        if(iblks==1)nremesh=0
        call TIME(char_time)
        print *, 'time: ', char_time
        write(chkunit,*)'time: ', char_time

        write(7,*)'ninit=',ninit,'iblks=',iblks,'restart=',restart
        if(ninit/=0.and.iblks==1.and.restart==0)call read_initial

        if (uinitial(iblks)==1) then
            result_zero=0.0
            if(allocated(result_first))result_first=0.0
            if(allocated(result_second))result_second=0.0
        endif
        print *,'a1'

        appear_p=appear
        do igroup=1,ngroup
            appear(igroup)=appear_process(igroup,iblks)
            group(igroup)%matno=matno_process(igroup,iblks)
            if(appear_process(igroup,iblks)==0.and.              &
                appear_process(igroup,iblks-1)==1)                &
                appear(igroup)=-1
        end do

        !20231215YL
        if(restart==0) then !20231008
            if(gamamax/=0) then !yuanli20230926 等效线性化土体动力本构
                open(gamamaxunit,file=probn(1:len1)//'.gamax')
                call readgamamax
            endif
        else
            gamamax=0
        endif

        if(ninistn>0.and.stnunit/=0)call read_permanent_strain !20231010
        !20231215YL

        ! transform the value in the old mesh to the new mesh
        ! meshc---1, coarse mesh to fine mesh, for stress concentration problem
        !         2, fine mesh to coarse mesh, for temperature and creep problem
        !         3, mixed, for fluid dynamic problem with consideration of free surface

        !! the following is to modify the interpolation matrix
        if (meshc==1.or.meshc==2)then  !2003/10/31
            appear=appear_process(:,iblks)
            if (meshc==2)then
                do igroup=ngroup0+1,ngroup
                    cgroup=group(igroup)%cgroup
                    appear(-cgroup)=0
                end do
                goto 100
            endif
            !!!!!! obtain the average result to set to the refined or coarse mesh
            do igroup=1,ngroup0
                cgroup=abs(group(igroup)%cgroup)
                if(cgroup==0) goto 10
                if (appear(igroup)==1)then
                    call judge_fine_mesh(igroup,icjr)
                    if (icjr==1)then
                        appear(igroup)=0
                        appear_process(igroup,iblks:nblks)=0
                        appear(cgroup)=1
                        appear_process(cgroup,iblks:nblks)=1
                    endif   ! for conditions of cgroup
                endif
10              continue
            end do ! for igroup
        endif ! for meshc

100     if (meshc==1.or.meshc==2)then  ! 2003/10/31
            write(7,992)appear(1:ngroup0)
            write(7,992)appear(ngroup0+1:ngroup)
992         format(20i5)
            do itotv=1,ntotv
                if(trans(itotv)%nintf/=0)deallocate(trans(itotv)%rintf,trans(itotv)%listf)
                trans(itotv)%nintf=0
            end do
            allocate(dinterp(npoin))
            dinterp=0
            allocate(icpoi(npoin))
            if (meshc==2) then
                icpoi=0
                do igroup=1,ngroup0
                    cgroup=abs(group(igroup)%cgroup)
                    if (cgroup==0.and.appear_process(igroup,iblks)==1)then
                        DO ielgroup = 1,group(igroup)%nelgroup
                            ielem = group(igroup)%list(ielgroup)
                            lnods =>element(ielem)%field(1)%lnods_f
                            icpoi(lnods)=1
                            nullify(lnods)
                        end do
                    endif
                end do

                do igroup=ngroup0+1,ngroup
                    cgroup=abs(group(igroup)%cgroup)
                    if (appear_process(igroup,iblks)==1.and.appear_process(cgroup,iblks)==1)then
                        DO ielgroup = 1,group(igroup)%nelgroup
                            ielem = group(igroup)%list(ielgroup)
                            lnods =>element(ielem)%field(1)%lnods_f
                            icpoi(lnods)=1
                            nullify(lnods)
                        end do
                    endif
                end do
                do ipoin=1,npoin
                    nintf=trans_c(ipoin)%nintf
                    if (icpoi(ipoin)==1.or.nintf<0)then
                        dinterp(ipoin)=1
                    endif
                end do
            elseif(meshc==1)then
                do igroup=ngroup0+1,ngroup
                    if (appear_process(igroup,iblks)==1)then
                        icpoi=0
                        DO ielgroup = 1,group(igroup)%nelgroup
                            ielem = group(igroup)%list(ielgroup)
                            lnods =>element(ielem)%field(1)%lnods_f
                            icpoi(lnods)=1
                            nullify(lnods)
                        end do
                        dinterp=dinterp+icpoi
                    endif
                end do
            endif

            deallocate(icpoi)

            do ipoin=1,npoin

                nintf=abs(trans_c(ipoin)%nintf)
                if ((meshc==2.and.nintf/=0.and.dinterp(ipoin)==1).or.   &
                    (meshc==1.and.trans_c(ipoin)%nintf>0.and.dinterp(ipoin)==1))then  !x2

                    !write(chkunit,*)'ipoin=',ipoin,'nintf=',nintf
                    !write(chkunit,*)'list=',trans_c(ipoin)%listf(:)
                    do idofn=1,cdofn   !4
                        itotv=nodfn(idofn,ipoin)
                        if (itotv/=0) then  !x3
                            trans(itotv)%nintf=nintf
                            allocate(trans(itotv)%listf(nintf),trans(itotv)%rintf(nintf))
                            if(meshc==1)result_zero(itotv)=0.
                            if (meshc==2.or.meshc==1)then
                                if(allocated(result_first))  result_first(itotv)=0.0
                                if(allocated(result_second))result_second(itotv)=0.0
                                if(allocated(torel))                torel(itotv)=0.0
                                if(allocated(toforl))              toforl(itotv)=0.0
                            endif

                            do jnode=1,nintf  !3
                                jpoin=trans_c(ipoin)%listf(jnode)
                                dfact=trans_c(ipoin)%rintf(jnode)
                                jtotv=nodfn(idofn,jpoin)
                                trans(itotv)%listf(jnode)=jtotv
                                trans(itotv)%rintf(jnode)=dfact
                                if(meshc==1)result_zero(itotv)=result_zero(itotv)+result_zero(jtotv)*dfact
                                if (meshc==1.or.meshc==2)then
                                    if(allocated(result_first))  result_first(itotv)=result_first(itotv)+dfact*result_first(jtotv)
                                    if(allocated(result_second))result_second(itotv)=result_second(itotv)+dfact*result_second(jtotv)
                                    if(allocated(torel)) torel(itotv) =torel(itotv) +dfact*torel(jtotv)
                                    if(allocated(toforl))toforl(itotv)=toforl(itotv)+dfact*toforl(itotv)
                                endif
                            enddo !3
                        endif   !x3
                    end do  !4
                endif     !x2
            end do
            deallocate(dinterp)
        endif

        write(7,*)'ikindks=',ikindks
        if(ikindks/=0)call steel_spring_parameter  !!steel 2006


        call prescrib_set  !20221124

        if(nbackdT==1)then !20210820

            do igapbf=1,nbackf
                igapb=backf(igapbf)%groupb
                npgblock=gapb(igapb)%npgblock
                do i0=1,npgblock
                    igaps=gapb(igapb)%nodegblock_igaps(i0)
                    ipairs=gapb(igapb)%nodegblock_ipairs(i0)
                    ij=gapb(igapb)%nodegblock_onetwo(i0)
                    ipoin=gaps(igaps)%pairnode(ij,ipairs)
                    do idofn=1,mdofn
                        if(lmdofn(idofn)==0)cycle  !20230430
                        itotv=nodfn(lmdofn(idofn),ipoin)
                        if(itotv==0)cycle
                        iffix(itotv)=6  !20210820
                        fixed(itotv)=0. !20210820
                    end do
                end do
            end do
        endif

        call external_load_2
        if(block_stab==0) &   !20200331
            call contact_pair_process !ctt2005  !zhao 2007.04.16
        call boundt !! temperature
        if(ADINA/=0.and.iblks==runblks)call GHM2ADINA
        if(ADINA==1.and.iblks<runblks)cycle
        if(ADINA==1.and.iblks==runblks)stop 'stop for ADINA==1!'
        line_load_block(iblks)=lineload
        line_temp_block(iblks)=linet  !! temperature
        call TIME(char_time)
        print *, 'time: ', char_time
        write(chkunit,*)'time(solve): ', char_time
        operation='SET'
        call solve
        call TIME(char_time)
        print *, 'time: ', char_time
        write(chkunit,*)'time(solve_set): ', char_time
        if(nlayer==2.and.type_solver=='PROFILE'.and.solver_iter==2)call nonzero_stiff_pcg
        call TIME(char_time)
        print *, 'time: ', char_time
        write(chkunit,*)'time(solve_pcg): ', char_time
        call modf_time_order ! zhao 2007/04/10
        if (type_problem=='Q') then

            if(relis==1.and.block_stab==0)then
                call STATIC_U_reli
            elseif(relis==1.and.block_stab>=1)then
                call static_rigid_reli
            elseif(block_stab>=1.and.ebody==0)then
                call static_rigid_1
            else if(mdofn==7)  then
                if(lmdofn(7)/=0)call static_U_P
            else if(mdofn==8)  then
                print *,' to static_u_pw'
                if(lmdofn(8)/=0)call static_U_Pw
            else if(nbackf/=0)then
                if(nbackdT==0)then
                    call back_analysis  !20150925
                else
                    call back_d_analysis  !20210820
                endif
            else
                call static_U
            endif
        endif
        !if(type_problem=='Q'.and.mdofn>=10.and.lmdofn(10)/=0)call modf_time_order  !zhao 05/07/22
        !call modf_time_order ! zhao 2007/04/10
        if (type_problem/='Q'.and.type_problem/='E'.and.type_problem/='W') then !freq2006
            call modf_time_order
            if(type_solver=='EXPLICIT') then
                call explicit
            else
                call time_dependent
            endif
        elseif(type_problem=='E') then
            call response_spectrum
        endif

        if(type_problem=='W') call frequency_analysis !freq2006

        call gpvar1_initial                                         !! 4/5/98

        restart=0
        lincs=0
        if(meshc/=0) appear_p=appear  !2003/10/31

        if (rmesh/=0)then  !!2004/7/12

            cmesh=1
            if (type_load=='DISCONTROL')then
                do itcurve=1,ntcurve
                    type_curve=tcurves(itcurve)%type_curve
                    if (type_curve=='DISCONTROL')then
                        ic=tcurves(itcurve)%ttime_curve(1)
                        yload=abs(prescrib(ic)%rdofix)
                        err=(yload-xload)/yload
                        if(abs(err)<.02.or.err<-.05) cmesh=0
                    endif
                end do
            end if

            print *,'cmesh=',cmesh
            nremesh=nremesh+1
            print *,'nremesh=',nremesh
            if(nremesh>1)call result_store_of_fine_mesh
            if (cmesh==1)then
                nelc=0;nelc1=0
                call find_remesh_stran
                !call find_remesh_element1
                !if(rmesh==2) &
                !call find_remesh_element2
                !write(7,*)'nelc=',nelc,'nelc1=',nelc1
                !if(nelc/=0) &
                !write(7,*)'listnelc=',listnelc
                call mesh_refine(nremesh)
                if(nelc/=0)deallocate(listnelc)
                if(nelc1/=0)deallocate(listnelc1)
            endif
        endif !!2004/7/12
    end do    !iblks

    !20231215YL
    !&  等效线性动剪应力输出 yuanli20230926 !20231008
    if (restart==0.and.gamamax/=0) then
        call writegamamax
    endif

    if(type_problem=='F')then !20231009
        call out_gid_max
    endif
    !20231215YL


    end subroutine process_analysis  !20190810

    SUBROUTINE stab_initialize

    character(10) fieldid,class,name,material
    integer(ink) igroup,ielgroup,ielem,matno,order_int,index,icr

    if(ninit/=0)stop 'stop for ninit/=0!!'

    result_zero=0.

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
                icr=0
                if(material=='CONCRETE')icr=props(matno)%mechanical%solid%Concrete%icr
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    if(ice0(ielem)==1) goto 100
                    element(ielem)%field(1)%gpvar0=0.
                    element(ielem)%field(1)%gpvar=0.
                    !!contact
100                 continue
                end do
            endif
        endif !--2
    end do !--1

    end SUBROUTINE stab_initialize
