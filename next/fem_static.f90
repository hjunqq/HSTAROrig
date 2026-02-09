    SUBROUTINE STATIC_U

    logical logx
    character(80)text
    integer(ink) i,itotv,ielem,irst,trstep0,ipoin,idofn,ij,ij0,idofix,ldofix,idelgroup,i0,ipairs
    real   (irk) xtime,time_begin,detal,ttime0,coef
    real   (irk),allocatable::rvectorm(:),value(:)
    integer(ink) iintf,nintf,iieq,njntf,bblks   !!int2000
    integer(ink) iincs_i,iblks_i,istep_i,inode,jnode,ivalue,ivalue_point  !20200819
    integer(ink),pointer::listf(:)  !20200819
    real   (irk),pointer::rintf(:)  !20200819

    integer(ink) igapb,npgblock,jpoin,igaps,ipair,idimn,itotvbt,jdimn, &  !!ctt2005
        jtotv,kpoin,lpoin,jtotvbt,npairs,cwater,jgaps,jpair,kdimn   !!ctt2005
    real   (irk),allocatable::rot(:,:),tofor0(:)  !!ctt2005
    real   (irk),allocatable::unitl(:),unitg(:),cmatrixl(:,:) !!ctt2005

    if (meshc==1.or.rmesh/=0)rewind(mainunit)
    if(Bparameter/=0.and.iblks==1)rewind(mainunit)  !20230902
    if(Bparameter/=0.and.iblks==1)rewind(upliftunit)


    read(mainunit,*)text
    read(mainunit,*)nincs

    print *,' in static_U**'

    if(ngaps/=0.or.nrcsteel/=0)allocate(tofor0(ntotv)) !!ctt2005

    do iincs=1,lincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta,cwater,Qstatic
        read(mainunit,*)toler_force,toler_var(1:mdofn)
        if(cwater/=0.and.delgroup>0)then
            do idelgroup=1,delgroup
                read(mainunit,*)text
            end do
        end if
        if(Qstatic/=0) then !20221104
            do i0=1,6
                read(mainunit,*)text
            end do
        endif

    end do

    xtime=0.0
    !write(7,*)'iffix(nodfn(1:3,25))=',iffix(nodfn(1:3,25))

    do iincs=lincs+1,nincs
        print *,'iincs=',iincs,'lincs=',lincs

        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta,cwater,Qstatic
        read(mainunit,*)toler_force,toler_var(1:mdofn)
        if(cwater/=0.and.delgroup>0)then
            allocate(coef_water(delgroup,nstep))
            do idelgroup=1,delgroup
                read(mainunit,*)i0,coef_water(idelgroup,:)
            end do
        end if


        if(Qstatic/=0) then !20221104
            read(mainunit,*)text  !20221104
            allocate(qstatic_force)  !20221104
            allocate(qstatic_force%appearg(ngroup),qstatic_force%qfactor(ndimn),qstatic_force%cor_coef(2,Qstatic))
            read(mainunit,*)qstatic_force%iaxe
            read(mainunit,*)qstatic_force%appearg
            read(mainunit,*)qstatic_force%qfactor
            read(mainunit,*)qstatic_force%cor_coef(1,:)
            read(mainunit,*)qstatic_force%cor_coef(2,:)
        endif !20221104

        print *,'cwater,Qstatic=',cwater,Qstatic

        ttime0=ttime
        trstep0=trstep
        do istep=inc_step,nstep,inc_step

            if(iblks>=stab_matde)call stab_initialize

            write(chkunit,*)'Increment step=',istep
            print *,'istep=',istep
            if(outintr>0.and.iblks>=outintr)trstep=trstep0+istep !20200226
            xtime=ditime*istep
            ttime=ttime0+ditime*istep !! only for output

            call dfact_time_curve(ttime)
            call modf_var_prescribed

            call saturation_judge  !20220409

            print *,'tcurvegravity=',tcurvegravity
            call gravity
            if(rmesh>0)call gravity1
            if(rmesh>1)call gravity2
            write(7,*)'cwater=',cwater,'delgroup=',delgroup
            if(cwater/=0.and.delgroup/=0)call step_water_pressure  !2013/3/18

222         call force_external


            !if(iblks/=1)mdiv=1   !5
            if(type_load=='LOAD2')mdiv=2  !!806
            do idiv=1,mdiv
                !! temperature
                if(type_load=='LOAD2'.and.idiv==2) goto 71
                call load_of_creep_and_temperature
                call creep_strain_of_rock_fill    !20130510
                call wetting_strain_of_rock_fill  !20220409
71              if(mdiv/=1)toform=toforl+(tofor-toforl)*idiv/mdiv
                if(type_load=='DISCONTROL')preact0=prescrib(1)%rdofix
                if((ngaps/=0.or.nrcsteel/=0).and.mdiv==1)tofor0=tofor !!ctt2005
                if((ngaps/=0.or.nrcsteel/=0).and.mdiv/=1)tofor0=toform !!ctt2005
                deltafi=0.0
                do igapb=1,ngapb !fzx  tcl
                    if(gapb(igapb)%nrdof==0)cycle
                    gapb(igapb)%rdisp_deltafi=0.
                enddo

                if(submodel==1.and.idiv==1)call value_submodel_boundary  !20210321


                do iiter=1,miter
                    iccontact=0 !zhao 05/07/30
                    !print *,'iblks=',iblks,'idiv=',idiv,'iiter=',iiter
                    if (istep==inc_step.and.iiter==1)then
                        call local_stress
                        call contact_state(0)
                        !call crack_state !crack 2006
                    endif

                    call algort

                    if (iiter==1.or.(kstat==2.and.iiter.le.2))then
                        !deltafi_ssorpbcg=deltafi !ssorpbcg
                        !deltafi=0.0
                        delitfi=0.0
                        call predict
                        do ielem=1,nelem   !!simo_rifai
                            if(associated(element(ielem)%alfa))element(ielem)%alfa=0.
                        end do  !!simo_rifai
                    endif


                    if(ikindks/=0) call strain_for_steel_bar !steel 2008

                    call stran0_creep4   !20180630  (博格斯模型蠕变初应变增量，因为应力增量在变化，所以每一迭代步求解，只适用于NSOLN=5）
                    if(iiter==1)   call effect_stres_modul_for_steel_beam !20211125
                    if(iiter==1)   call stiffness_for_bolt_spring  !20211125

                    if (nlayer/=2)then

                        if (kresl/=0.or.kthmat/=0) then
                            if(kresl/=0)call stiff_u
                            if(kresl/=0.and.rmesh>0.and.nelem1>0)call stiff_u1
                            if(kresl/=0.and.rmesh>1.and.nelem2>0)call stiff_u2
                            if(neuman==1.and.((kstat==2.and.iiter==2).or.&
                                (kstat/=2.and.istep==inc_step.and.iiter==1)))call write_stiff_u
                            if(kthmat/=0)call htmatrx

                            if(type_solver=='PROFILE'.and.   &
                                (neuman==1.and.((kstat/=2.and.(istep/=1.or.iiter/=1)).or.(kstat==2.and.iiter.gt.2))))goto 1
                            if(type_solver/='JPCG')global_stiff1=0.0
                            if(nonsym/=0.and.type_solver=='PROFILE')global_stiff2=0.0
                            if (type_solver=='JPCG'.and.outintr==0) then
                                do ielem=1,nelem
                                    element(ielem)%estif=0.0
                                end do
                            endif
                            call estif_assemble
                            if(nbspring>0) &       !20150925
                                call assemble_back_spring  !20150925
                            if(ground_inf/=0) call semi_inf_space_assemble  !20231010

                            if(nonsym==0)then !20240312 YL
                                do itotv=1,ntotv
                                    if (totveq(itotv)/=0)then
                                        if(abs(global_stiff1(iseq(totveq(itotv)))).le.1.e-5)global_stiff1(iseq(totveq(itotv)))=1.e30
                                    endif
                                enddo
                            endif !20240312 YL

                        endif

                    else !if (nlayer/=2)then

                        if(kresl_layer1/=0.or.kresl_layer2/=0)call stiff_u
                        print *,'kresl_layer=',kresl_layer1,kresl_layer2
                        if(kresl_layer1/=0)global_stiff1(1:iseq(neq_layer1))=0.
                        if(kresl_layer2/=0)global_stiff1(iseq(neq_layer1)+1:iseq(neq))=0.
                        if(nonsym==1.and.kresl_layer1/=0)global_stiff2(1:iseq(neq_layer1))=0.
                        if(nonsym==2.and.kresl_layer2/=0)global_stiff2(iseq(neq_layer1)+1:iseq(neq))=0.

                        call estif_assemble

                        if(nonsym==0)then !20240312 YL
                            do itotv=1,ntotv
                                if (totveq(itotv)/=0)then
                                    if(abs(global_stiff1(iseq(totveq(itotv)))).le.1.e-5)global_stiff1(iseq(totveq(itotv)))=1.e30
                                endif
                            enddo
                        endif !20240312 YL

                    endif !if (nlayer/=2)

1                   if(type_load=='LOAD2'.or.(kstat==2.and.iiter.le.2).or.(type_load/='LOAD2'.and.kstat/=2.and.iiter==1).or.  &
                        (ngaps/=0.and.istatec==0)) then	  !! for temperature 20130510
                        !if(type_load=='LOAD2'.and.idiv==2)then  !20130510  !20220607
                        ! deltafi=0.0
                        ! delitfi=0.0
                        ! endif
                        call gpvar2_initial
                        if (ninit/=0.and.(kinit==2.and.iincs==1.and.istep==1)) then  !20201203
                            call eload_initialize
                            call eload_initial_stress
                            call force_release
                            where(totveq==0)
                                torel=0.0
                            endwhere
                        endif !20201203


                        call eload_initialize
                        write(7,*)'residu_f1'
                        call residu_f  !!!!!!

                        if(rmesh>0.and.nelem1>0)call residu_f1
                        if(rmesh>1.and.nelem2>0)call residu_f2
                        call eload_field
                        if(nbspring>0) & !20150925
                            call eload_back_spring  !20150925
                        if(ground_inf/=0)call semi_inf_load
                        call force_internal
                    endif !for iiter==1 and istep==inc_step .and.idiv==1  temperature
                    if (mdiv/=1) then
                        if(idiv==1.and.iiter==1.and.allocated(torel))toform=toform+torel
                    else
                        if(iiter==1.and.allocated(torel))tofor=tofor+torel
                    endif

                    if(ngaps/=0.and.iblks>=iblks_bt.and.mdiv==1)call ctfor_to_tofor(tofor0,tofor)  !!ctt2005
                    if(ngaps/=0.and.iblks>=iblks_bt.and.mdiv/=1)call ctfor_to_tofor(tofor0,toform)  !!ctt2005

                    if(nrcsteel/=0.and.mdiv==1)call csfor_to_tofor(tofor0,tofor)  !!20210328
                    if(nrcsteel/=0.and.mdiv/=1)call csfor_to_tofor(tofor0,toform)  !!20210328


                    if(neuman==1.and.(istep/=1.or.iiter/=1).and.kresl/=0)  goto 2  !ctt2005 , change position!
                    if(type_nl==8.and.(iiter>1.or.(kstat==2.and.iiter>2))) goto 2  !MNR
                    if ((type_solver=='PROFILE'.or.type_solver=='PARDISO').and.kresl/=0)then
                        operation='FACTORIZE'
                        call solve
                    end if
2                   continue

                    logx=ngaps/=0.and.(iiter==1.and.istep==inc_step).and.iincs==(lincs+1) !20200331
                    if (logx)then !ctt2005
                        if (restart_ctt==0)then !restart_ctt
                            kdimn=ndimn
                            if(block_stab==1)kdimn=3*(ndimn-1) !2015/11/17
                            allocate(rot(kdimn,kdimn))
                            rot=0.
                            do igapb=1,ngapb
                                if(block_appear_process(igapb,iblks)==0)cycle    !20200331
                                npgblock=gapb(igapb)%npgblock
                                gapb(igapb)%cmatrix=0.

                                if(gapb(igapb)%eblock==0) cycle  !2017/11/19
                                do ipoin=1,npgblock
                                    igaps=gapb(igapb)%nodegblock_igaps(ipoin)
                                    ipair=gapb(igapb)%nodegblock_ipairs(ipoin)
                                    ij=gapb(igapb)%nodegblock_onetwo(ipoin)


                                    coef=1.
                                    if(ij==2)coef=-1.
                                    rot(1:ndimn,1:ndimn)=gaps(igaps)%rot(:,:,ipair)
                                    if(kdimn>ndimn)then
                                        if(ndimn==2)rot(3,3)=1.
                                        if(ndimn==3)rot(4:6,4:6)= rot(1:ndimn,1:ndimn)
                                    endif

                                    allocate(unitl(kdimn),unitg(kdimn))
                                    do idimn=1,kdimn

                                        itotvbt=(ipoin-1)*kdimn+idimn
                                        unitl=0.
                                        unitl(idimn)=1.*coef


                                        unitg=transpose(rot).x.unitl
                                        rvector=0.
                                        call  unit_force_trans(igapb,kdimn,ij,unitg,igaps,ipair,rvector)

                                        operation='SOLVE'
                                        call solve


                                        do kpoin=1,npgblock
                                            jgaps=gapb(igapb)%nodegblock_igaps(kpoin)
                                            jpair=gapb(igapb)%nodegblock_ipairs(kpoin)
                                            ij0=gapb(igapb)%nodegblock_onetwo(kpoin)         !2017/04/03
                                            call result_node_to_center(kdimn,ij0,jgaps,jpair,result,unitg)

                                            do jdimn=1,kdimn
                                                jtotvbt=(kpoin-1)*kdimn+jdimn
                                                gapb(igapb)%cmatrix(jtotvbt,itotvbt)=unitg(jdimn)
                                            end do
                                        end do  !kpoin
                                    end do  !idimn
                                    deallocate(unitl,unitg)
                                end do  !ipoin

                                do ipoin=1,npgblock
                                    igaps=gapb(igapb)%nodegblock_igaps(ipoin)
                                    ipair=gapb(igapb)%nodegblock_ipairs(ipoin)
                                    ij=gapb(igapb)%nodegblock_onetwo(ipoin)
                                    coef=1.
                                    if(ij==2)coef=-1.
                                    rot=0.
                                    rot(1:ndimn,1:ndimn)=gaps(igaps)%rot(:,:,ipair)

                                    if(kdimn>ndimn)then
                                        if(ndimn==2)rot(3,3)=1.
                                        !if(ndimn==3)rot(1:3,4:6)= rot(1:ndimn,1:ndimn)
                                        if(ndimn==3)rot(4:6,4:6)= rot(1:ndimn,1:ndimn)
                                        !if(ndimn==3)rot(4:6,1:3)= rot(1:ndimn,1:ndimn)
                                    endif

                                    rot=coef*rot

                                    allocate(cmatrixl(kdimn,kdimn))
                                    do jpoin=1,npgblock
                                        jtotv=(jpoin-1)*kdimn
                                        itotv=(ipoin-1)*kdimn
                                        cmatrixl=gapb(igapb)%cmatrix(itotv+1:itotv+kdimn,jtotv+1:jtotv+kdimn)
                                        gapb(igapb)%cmatrix(itotv+1:itotv+kdimn,jtotv+1:jtotv+kdimn)=   &
                                            rot.x.cmatrixl
                                    end do
                                    deallocate(cmatrixl)
                                end do
                            end do  !igapb
                            call forAdirect !fzx !形成A矩阵

                            do igapb=1,ngapb
                                !write(7,*)'igapb=',igapb,'ntotv_bt=',gapb(igapb)%ntotv_bt,'camatrix='
                                do itotvbt=1,gapb(igapb)%ntotv_bt
                                    !write(7,*)gapb(igapb)%cmatrix(itotvbt,:)
                                    do jtotvbt=1,gapb(igapb)%ntotv_bt
                                        write(recttunit)gapb(igapb)%cmatrix(itotvbt,jtotvbt)
                                    enddo
                                enddo
                            enddo  !igapb

                            deallocate(rot)

                        elseif(restart_ctt==1)then !restart_ctt
                            call forAdirect !fzx !形成A矩阵
                            rewind(recttunit)
                            do igapb=1,ngapb
                                npgblock=gapb(igapb)%npgblock
                                do itotvbt=1,gapb(igapb)%ntotv_bt
                                    do jtotvbt=1,gapb(igapb)%ntotv_bt
                                        read(recttunit)gapb(igapb)%cmatrix(itotvbt,jtotvbt)
                                    enddo
                                enddo
                            enddo
                        else !restart_ctt
                            write(*,*)'no such restart_ctt!!'
                            stop
                        endif !restart_ctt
                    endif  !!ctt2005
                    !logx=nrcsteel/=0.and.(iiter==1.and.istep==inc_step).and.iincs==(lincs+1) !20220311
                    !if (logx) call cmatrix_c_formation !20210328

                    if(kresl/=0) & !20220311
                        call cmatrix_c_formation !20220311



90                  format(10e12.5)
                    rvector=0.0
                    !write(chkunit,*)'iiter=',iiter,'itotv ieq  tofor stfor'
                    !write(chkunit,*)'itotv=','ieq=','stif=','rvector='

                    if (type_solver/='JPCG') then
                        do itotv=1,ntotv
                            !if(abs(tofor(itotv)-stfor(itotv))>1.e-5) &

                            if (totveq(itotv)/=0)then
                                if (mdiv/=1)then
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                        toform(itotv)-stfor(itotv)

                                else
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                        tofor(itotv)-stfor(itotv)
                                    !if(abs(rvector(totveq(itotv)))>1.e-8)write(7,*)itotv,tofor(itotv),stfor(itotv)
                                    !write(7,*)itotv,tofor(itotv),stfor(itotv)
                                endif
                                if(nonsym==0)then !20240312 YL
                                    if(abs(global_stiff1(iseq(totveq(itotv)))).le.1.e-5)global_stiff1(iseq(totveq(itotv)))=1.e20
                                endif !20240312 YL
                                !if(abs(rvector(totveq(itotv)))>1.e-8) &
                                !write(chkunit,52)itotv,totveq(itotv),global_stiff1(iseq(totveq(itotv))),rvector(totveq(itotv))


                            endif
                            !write(7,*)itotv,tofor(itotv),stfor(itotv)
                        end do
                        !stop
                        !!int2000

                        do itotv=1,ntotv
                            nintf=trans(itotv)%nintf
                            if (nintf/=0) then
                                do iintf=1,nintf
                                    iieq=totveq(trans(itotv)%listf(iintf))
                                    if(iieq/=0)rvector(iieq)=rvector(iieq)+  &
                                        (tofor(itotv)-stfor(itotv))*trans(itotv)%rintf(iintf)
                                end do
                            endif
                        end do

                        !!int2000

                    else !if (type_solver/='JPCG') then

                        if(mdiv/=1)rvector=toform-stfor
                        if(mdiv==1)rvector=tofor -stfor
                    endif
52                  format(2I10,3e15.5)



                    if (type_load=='ARCLENGTH'.and.kresl/=0) then
                        allocate(rvectorm(neq))
                        rvectorm=rvector
                        rvector=0.0
                        if (type_solver/='JPCG') then
                            do itotv=1,ntotv
                                if(totveq(itotv)/=0) &
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+tofor_arclength(itotv)
                            end do
                        else
                            rvector=tofor_arclength
                        endif
                        operation='SOLVE'
                        call solve
                        delta_arclength=result
                        rvector=rvectorm
                        deallocate(rvectorm)
                    endif

                    !write(7,*)'idiv=',idiv,'iiter=',iiter
                    !write(7,*)'rvector=',rvector

                    if (type_nl==8)then
                        if(kstat/=2)call bfgsr(iiter)
                        if(kstat==2)call bfgsr(iiter-1)
                    else
                        operation='SOLVE'
                        call solve
                    endif

                    !write(7,*)'result=',result

                    if(ngaps/=0.and.iblks>=abs(iblks_bt))call solve_ctt  !!ctt2005
                    if(nrcsteel/=0) call solve_bond_force_of_cs !20210328

                    !
                    !write(7,*)'tofor***'
                    !do itotv=1,ntotv
                    !    if(abs(tofor(itotv))>1.e-3) &
                    !    write(7,*)itotv,tofor(itotv)
                    ! end do

                    if (type_load=='ARCLENGTH')then
                        call find_dfact_of_arclength (irst)
                        if (irst==1) then
                            time_begin=tcurves(arc_curve)%time_begin
                            detal=tcurves(arc_curve)%detal
                            if (abs(ttime-time_begin-ditime).le.1.e-8)detal=tcurves(arc_curve)%fact_inc
                            tcurves(arc_curve)%detal=detal*.5
                            if (abs(ttime-time_begin-ditime).le.1.e-8)tcurves(arc_curve)%fact_inc=detal*.5
                            goto 222
                        endif
                    endif

                    if(neuman==1.and.(istep/=1.or.iiter/=1).and.kresl/=0)call neuman_expan

                    call TIME(char_time)
                    print *, 'time: ', char_time
                    write(chkunit,*)'time: ', char_time


                    !write(7,*)'varupdate'
                    call varupdate
                    !write(7,*)'af varupdate'
                    call relative_dis_watertight !20231007 止水 !20240305
                    call eload_initialize

                    if(ikindks/=0) call strain_for_steel_bar !steel 2008
                    !write(7,*)'bbxx residu_f'
                    write(7,*)'residu_f2'
                    call residu_f

                    !write(7,*)'aaxx residu_f'
                    if(rmesh>0.and.nelem1>0)call residu_f1
                    if(rmesh>1.and.nelem2>0)call residu_f2
                    if(type_load/='LOAD2'.or.(type_load=='LOAD2'.and.idiv==2))then   !806
                        call eload_field
                        if(nbspring>0) & !20150925
                            call   eload_back_spring  !20150925
                        if(ground_inf/=0)call semi_inf_load

                        call reaction_prescribed

                        call conver_load
                        if(nchek==0) call conver_nodal_value
                        if(type_load=='ARCLENGTH')tcurves(arc_curve)%piter=iiter

                        !!!!!!!!!!!!!!!!
                        !call local_stress
                        !call contact_state(1)

                        !**************************

                        !if(nchek==0.and.iccontact==1)exit
                        !!!!!!!!!!!!!!!
                        if(nchek==0)exit !tcl
                    endif !806
                    if(type_load=='LOAD2'.and.idiv==1) goto 10
                    print *,'miter=',miter,'iiter=',iiter
                end do   !! loop for iiter
10              continue

                if(type_load/='LOAD2'.or.(type_load=='LOAD2'.and.idiv==2))then   !20220607
                    call local_stress
                    call contact_state(1)
                    !write(7,*)'gpvar1=',element(1)%field(1)%gpvar(1:6,1)
                    !if(type_load/='LOAD2') then   !20220607
                    if(istatec==0) &
                        call state_and_stiff_2021
                    !write(7,*)'icttstif_static_u=',icttstif
                    do igaps=1,ngaps
                        npairs=gaps(igaps)%npairs
                        do ipairs=1,npairs
                            if(gaps(igaps)%pair_process(ipairs)==0)cycle  !20200331
                            gaps(igaps)%dxyz0(:,ipairs)=gaps(igaps)%dxyz(:,ipairs)
                            gaps(igaps)%ctforce0(:,ipairs)=gaps(igaps)%ctforce(:,ipairs)
                            gaps(igaps)%state0(ipairs)=gaps(igaps)%state(ipairs)
                            gaps(igaps)%damage0(ipairs)=gaps(igaps)%damage(ipairs)
                            gaps(igaps)%kxyz0(:,:,ipairs)=gaps(igaps)%kxyz(:,:,ipairs)
                        end do
                    end do
                    !endif   !20220607

                    !if(type_load/='LOAD2') &   !20220607
                    call gpvarupdate

                    if(rmesh>0.and.nelem1>0)call gpvarupdate1
                    if(rmesh>1.and.nelem2>0)call gpvarupdate2

                endif   !20220607
            end do    !! for idiv

            !if(type_load=='LOAD2') &   !20220607
            !call gpvarupdate   !20220607

            if(Blarge==1) then  !20221102
                call update_coord_blarge
                call modf_element_local_direction
            endif !20221102
            if(modf_dis_blocks(iblks)==1)call construction_dis_modify


            !if(type_load=='LOAD2') then   !20220607
            !if(istatec==0) &
            !              call state_and_stiff_2021
            !do igaps=1,ngaps
            !npairs=gaps(igaps)%npairs
            !do ipairs=1,npairs
            !     if(gaps(igaps)%pair_process(ipairs)==0)cycle  !20200331
            !    gaps(igaps)%state0(ipairs)=gaps(igaps)%state(ipairs)
            !   gaps(igaps)%damage0(ipairs)=gaps(igaps)%damage(ipairs)
            !   gaps(igaps)%ctforce0(:,ipairs)=gaps(igaps)%ctforce(:,ipairs)
            !   gaps(igaps)%dxyz0(:,ipairs)=gaps(igaps)%dxyz(:,ipairs)
            !   gaps(igaps)%kxyz0(:,:,ipairs)=gaps(igaps)%kxyz(:,:,ipairs)
            !end do
            !end do
            !    endif      !20220607


100         toforl=tofor

            if (istep/noutn*noutn==istep)then
                iwriten=iwriten+1
                call out_record
                call outputres !for output
            endif


            !if (kstab==0.) then
            if(nforce/=0.or.ngaps/=0)call force_interface
            !else
            if(kstab/=0.)call safety_factor
            !endif
            if (istep/noutf*noutf==istep)then
                !if(kstab==0.and.nforce/=0)call write_force_interface
                if(nforce/=0.or.ngaps/=0)call write_force_interface

                call out_full_write
                if(outplot(1:3)=='GID')   call OUT_GID_WRITE
                if(outplot(1:6)=='COSMOS')call OUT_COSMOS_WRITE
            endif

            if(istep/nresta*nresta==istep)call resta_read_write(-1)

            if(Bparameter>0)then !20230523
                !Value_observ(:)%value_computation=0.
                do ivalue=1,mvalue
                    !if(Value_observ(ivalue)%ic==0)cycle
                    iblks_i=Value_observ(ivalue)%iblks
                    iincs_i=Value_observ(ivalue)%iincs
                    istep_i=Value_observ(ivalue)%istep
                    idofn =lmdofn(Value_observ(ivalue)%idofn)
                    ivalue_point=Value_observ(ivalue)%ivalue_point
                    if(iblks_i==iblks.and.iincs_i==iincs.and.istep_i==istep)then
                        nintf=para_points(ivalue_point)%nintf
                        listf=>para_points(ivalue_point)%listf
                        rintf=>para_points(ivalue_point)%rintf
                        if(Bparameter==1)Value_observ(ivalue)%value_computation=dot_product(rintf,result_zero(nodfn(idofn,listf)))
                        !if(ivalue_point==169)then  !20230717
                        !    write(7,*)'ivalue_point=',ivalue_point,'ivalue=',ivalue
                        !    write(7,*)'result_zero=',result_zero(nodfn(idofn,listf))
                        !    write(7,*)'rintf=',rintf
                        !    write(7,*)'value_computation=',Value_observ(ivalue)%value_computation
                        ! endif

                        if(Bparameter==2)Value_observ(ivalue)%value_computation=dot_product(rintf,deltafi(listf))
                        nullify(listf,rintf)
                    endif
                end do
            endif   !20230523

            if(Bparameter<0)then !20200812
                tbstep=tbstep+1
                do i=1,nback_point
                    bblks=freedom_for_back(4,i)  !20230523
                    if(bblks>iblks)cycle !20230523

                    inode=freedom_for_back(1,i)
                    idofn=freedom_for_back(2,i)
                    jnode=freedom_for_back(3,i)
                    print *,'inode=',inode,'idofn=',idofn,'jnode=',jnode
                    itotv=nodfn(lmdofn(idofn),inode)
                    if(jnode/=0)jtotv=nodfn(lmdofn(idofn),jnode)
                    if(Bparameter==-1)then
                        Value_vc(i,tbstep,istoch)=result_zero(itotv)
                        if(jnode/=0)Value_vc(i,tbstep,istoch)=Value_vc(i,tbstep,istoch)-result_zero(jtotv)
                    elseif(Bparameter==-2)then
                        Value_vc(i,tbstep,istoch)=deltafi(itotv)
                        if(jnode/=0)Value_vc(i,tbstep,istoch)=Value_vc(i,tbstep,istoch)-deltafi(jtotv)
                    end if
                end do
            endif  !20200812
            if((bparameter>=1.and.bparameter<=2).and.balgor>=1) call dudx

        end do     !! loop for istep

        if(Qstatic/=0) then !20221104
            deallocate(qstatic_force%appearg,qstatic_force%qfactor,qstatic_force%cor_coef) !20221104
            deallocate(qstatic_force)  !20221104
        endif !20221104

        !if(kstab==0.and.(nforce/=0.or.ngaps/=0))call write_force_interface
        if(cwater/=0.and.delgroup>0)deallocate(coef_water)
    end do !!iincs
    !if(winit==-1) call out_next_write
    if(winit==-1*iblks) call out_next_write !20231215YULI

    !if(ngaps/=0)deallocate(tofor0)  !!ctt2005



    END SUBROUTINE STATIC_U

    SUBROUTINE STATIC_U_P

    character(80)text
    integer(ink) itotv,ielem,trstep0
    real   (irk) time,ttime0
    integer(ink) iintf,nintf,iieq   !!int2000

    print *,'in static_u_p'
    read(mainunit,*)text
    read(mainunit,*)nincs

    do iincs=1,lincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta
        read(mainunit,*)toler_force,toler_var(1:mdofn)
    end do

    time=0.0
    do iincs=lincs+1,nincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta
        read(mainunit,*)toler_force,toler_var(1:mdofn)
        ttime0=ttime
        trstep0=trstep
        do istep=inc_step,nstep,inc_step
            if (outintr>0.and.iblks>=outintr)trstep=trstep0+istep !20200226
            time=ditime*istep
            ttime=ttime0+ditime*istep        !! only for output
            call dfact_time_curve(ttime)
            call modf_var_prescribed
            call gravity
            call force_external
            do ielem=1,nelem   !!simo_rifai
                if (associated(element(ielem)%alfa))element(ielem)%alfa=0.
            end do  !!simo_rifai

            do iiter=1,miter

                if  (iiter==1) then
                    deltafi=0.0
                    delitfi=0.
                endif

                call algort

                if (iiter==1)call predict   ! new

                if (kresl/=0)call stiff_u
                if (ksmat/=0)call mcmatrx('P')
                if (kqmat/=0)call upwcouple

                if  (kresl/=0.or.ksmat/=0.or.kqmat/=0) then
                    if (type_solver/='JPCG')global_stiff1=0.0
                    if (nonsym/=0.and.type_solver=='PROFILE')global_stiff2=0.0
                    if  (type_solver=='JPCG')then
                        do ielem=1,nelem
                            element(ielem)%estif=0.0
                        end do
                    endif
                    call estif_assemble
                    call couple_assemble
                endif

                !            if (iiter==1.and.istep==inc_step) then   !! iiter==1 and istep==inc_step
                if (iiter==1) then   !! iiter==1 and istep==inc_step
1                   call eload_initialize
                    if ((ninit/=0.and.((iblks==1.and.kinit==1)).or.(kinit==2.and.iincs==1))) then

                        call eload_initial_stress
                        call gpvar2_initial
                        if (kinit==2.and.iincs==1)then
                            call force_release
                            where(totveq==0)
                                torel=0.0
                            endwhere
                        endif

                    endif
                    call residu_f
                    !           if(iblks/=1.or.(iblks==1.and.iincs/=1))call residu_f
                    call eload_field
                    call eload_couple

                    call force_internal

                endif    !! end for iiter==1 and istep==inc_step

                if (iiter==1.and.allocated(torel))tofor=tofor+torel


                rvector=0.0
                if (type_solver/='JPCG') then
                    do itotv=1,ntotv
                        if (totveq(itotv)/=0)     &
                            rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                            tofor(itotv)-stfor(itotv)
                    end do



                    !!int2000
                    do itotv=1,ntotv
                        nintf=trans(itotv)%nintf
                        if (nintf/=0) then
                            iieq=totveq(itotv)
                            if (iieq/=0)rvector(iieq)=0.
                            do iintf=1,nintf
                                iieq=totveq(trans(itotv)%listf(iintf))
                                if (iieq/=0) &
                                    rvector(iieq)=rvector(iieq)+  &
                                    (tofor(itotv)-stfor(itotv))*trans(itotv)%rintf(iintf)
                            end do
                        endif
                    end do
                    !!int2000
                else
                    rvector=tofor-stfor
                endif


                if (type_nl==8.and.(iiter>1.or.(kstat==2.and.iiter>2))) goto 2  !MNR

                if (type_solver=='PROFILE'.and.               &
                    (kresl/=0.or.ksmat/=0.or.khmat/=0.or.kqmat/=0)) then
                    operation='FACTORIZE'
                    call solve
                end if

2               if(type_nl==8) then
                    if (kstat/=2)call bfgsr(iiter)
                    if (kstat==2)call bfgsr(iiter-1)
                else
                    operation='SOLVE'
                    call solve
                endif

                where(iffix==0)
                    delitfi=result
                    deltafi=deltafi+delitfi
                endwhere
                where(iffix==0)
                    result_zero=result_zero+delitfi
                endwhere
                call eload_initialize
                call residu_f
                call eload_field
                call eload_couple
                call reaction_prescribed
                call conver_load

                if (nchek==0)call conver_nodal_value
                if (nchek==0)exit

            end do   !! loop for iiter

            call gpvarupdate

            if (istep/noutn*noutn==istep)then
                iwriten=iwriten+1
                call out_record
            endif

            if (istep/noutf*noutf==istep)then
                call out_full_write
                if (outplot(1:3)=='GID')   call OUT_GID_WRITE
                if (outplot(1:6)=='COSMOS')call OUT_COSMOS_WRITE
            endif
            if (istep/nresta*nresta==istep)call resta_read_write(-1)

        end do     !! loop for istep

    end do     !! loop for iincs

    END SUBROUTINE STATIC_U_P

    SUBROUTINE STATIC_U_PW


    logical logx
    character(80)text
    integer(ink) i,ii,itotv,ielem,irst,trstep0,ipoin,idofn,ij,ij0,idofix,ldofix,idelgroup,i0,ipairs
    real   (irk) xtime,time_begin,detal,ttime0,coef,pvalue,Tpredict !20230216
    real   (irk),allocatable::rvectorm(:),value(:)
    integer(ink) iintf,nintf,iieq,njntf,icdofn,ifixset   !!20230216  !!20220626
    integer(ink) iincs_i,iblks_i,istep_i,inode,jnode,ivalue,ivalue_point,bblks,mfixset,j  !20231130
    integer(ink),pointer::listf(:)  !20200819
    real   (irk),pointer::rintf(:)  !20200819

    integer(ink) igapb,npgblock,jpoin,igaps,ipair,idimn,itotvbt,jdimn, &  !!ctt2005
        jtotv,kpoin,lpoin,jtotvbt,npairs,cwater,jgaps,jpair,kdimn   !!ctt2005
    real   (irk),allocatable::rot(:,:),tofor0(:),midt(:)  !!ctt2005
    real   (irk),allocatable::unitl(:),unitg(:),cmatrixl(:,:) !!ctt2005
    real   (irk),allocatable::observstar(:),dtv(:),dtvi(:) !20230216
    real   (irk),allocatable::cmatrix_dtv(:,:),inv_cmatrix_dtv2(:,:)  !20230216



    print *,'static_U_Pw'
    if (meshc==1.or.rmesh/=0)rewind(mainunit)
    if(Bparameter/=0.and.iblks==1)rewind(mainunit)  !20230902
    if(Bparameter/=0.and.iblks==1)rewind(upliftunit)

    read(mainunit,*)text
    read(mainunit,*)nincs

    if(ngaps/=0.or.nrcsteel/=0)allocate(tofor0(ntotv)) !!ctt2005



    do iincs=1,lincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta,cwater,Qstatic
        read(mainunit,*)toler_force,toler_var(1:mdofn)
        if(cwater/=0.and.delgroup>0)then
            do idelgroup=1,delgroup
                read(mainunit,*)text
            end do
        end if
        if(Qstatic/=0) then !20221104
            do i0=1,6
                read(mainunit,*)text
            end do
        endif

    end do

    xtime=0.0
    do iincs=lincs+1,nincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta,cwater,Qstatic
        read(mainunit,*)toler_force,toler_var(1:mdofn)

        if(cwater/=0.and.delgroup>0)then
            allocate(coef_water(delgroup,nstep))
            do idelgroup=1,delgroup
                read(mainunit,*)i0,coef_water(idelgroup,:)
            end do
        end if


        if(Qstatic/=0) then !20221104
            read(mainunit,*)text  !20221104
            allocate(qstatic_force)  !20221104
            allocate(qstatic_force%appearg(ngroup),qstatic_force%qfactor(ndimn),qstatic_force%cor_coef(2,Qstatic))
            read(mainunit,*)qstatic_force%iaxe
            read(mainunit,*)qstatic_force%appearg
            read(mainunit,*)qstatic_force%qfactor
            read(mainunit,*)qstatic_force%cor_coef(1,:)
            read(mainunit,*)qstatic_force%cor_coef(2,:)
        endif !20221104


        ttime0=ttime
        trstep0=trstep

        do istep=inc_step,nstep,inc_step
            if(iblks>=stab_matde)call stab_initialize
            if (outintr>0.and.iblks>=outintr)trstep=trstep0+istep !20200226

            xtime=ditime*istep
            ttime=ttime0+ditime*istep        !! only for output

            call dfact_time_curve(ttime)
            call modf_var_prescribed

            call saturation_judge

            call gravity
            call loadfl

            !write(7,*)'cwater=',cwater,'delgroup=',delgroup
            if(cwater/=0.and.delgroup/=0)call step_water_pressure

222         call force_external

            if(type_load=='LOAD2')mdiv=2

            do idiv=1,mdiv

                if(type_load=='LOAD2'.and.idiv==2) goto 71
                call load_of_creep_and_temperature
                call creep_strain_of_rock_fill    !20130510
                call wetting_strain_of_rock_fill  !20220409
71              if(mdiv/=1)toform=toforl+(tofor-toforl)*idiv/mdiv
                if(type_load=='DISCONTROL')preact0=prescrib(1)%rdofix
                if((ngaps/=0.or.nrcsteel/=0).and.mdiv==1)tofor0=tofor !!ctt2005
                if((ngaps/=0.or.nrcsteel/=0).and.mdiv/=1)tofor0=toform !!ctt2005
                deltafi=0.0
                do igapb=1,ngapb !fzx  tcl
                    if(gapb(igapb)%nrdof==0)cycle
                    gapb(igapb)%rdisp_deltafi=0.
                enddo

                if(submodel==1.and.idiv==1)call value_submodel_boundary  !20210321
                do ielem=1,nelem   !!simo_rifai
                    if (associated(element(ielem)%alfa))element(ielem)%alfa=0.
                end do  !!simo_rifai

                do iiter=1,miter

                    iccontact=0
                    print *,'iblks=',iblks,'iincs=',iincs,'istep=',istep,'idiv=',idiv,'iiter=',iiter
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
                    if (istep==inc_step.and.iiter==1)then
                        call local_stress
                        call contact_state(0)
                        !call crack_state !crack 2006
                    endif
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

                    call algort
                    if (iiter==1.or.(kstat==2.and.iiter.le.2))then
                        delitfi=0.0
                        call predict
                        do ielem=1,nelem   !!simo_rifai
                            if(associated(element(ielem)%alfa))element(ielem)%alfa=0.
                        end do  !!simo_rifai
                    endif

                    if(ikindks/=0) call strain_for_steel_bar !steel 2008

                    call stran0_creep4   !20180630  (博格斯模型蠕变初应变增量，因为应力增量在变化，所以每一迭代步求解，只适用于NSOLN=5）
                    if(iiter==1)   call effect_stres_modul_for_steel_beam !20211125
                    if(iiter==1)   call stiffness_for_bolt_spring  !20211125

                    call porepr_w
                    if (kswkw/=0) call propty_w
                    if (kresl/=0             )call stiff_u
                    if (ksmat/=0.and.uwcpl==1)call mcmatrx('W')
                    if (kmass/=0)             call mcmatrx('U')
                    if (khmat/=0.and.uwcpl/=1)call hmatrx('W')
                    if (kqmat/=0.and.uwcpl/=0)call upwcouple
                    if (stabpw==1.and.khmat/=0)call stabpatch    !!stablize


                    if (kresl/=0.or.ksmat/=0.or.khmat/=0.or.kqmat/=0) then
                        if (type_solver/='JPCG')global_stiff1=0.0
                        if (nonsym/=0.and.type_solver=='PROFILE')global_stiff2=0.0
                        if (type_solver=='JPCG')then
                            do ielem=1,nelem
                                element(ielem)%estif=0.0
                            end do
                        endif
                        call estif_assemble
                        if (uwcpl==1)  call couple_assemble
                        if (stabpw==1) call stabpw_assemble           !! stablize
                    endif

                    if(nbspring>0) &       !20150925
                        call assemble_back_spring  !20150925
                    if(ground_inf/=0) call semi_inf_space_assemble

                    if(nonsym==0)then !20240312 YL
                        do itotv=1,ntotv   !20220618
                            if (totveq(itotv)/=0)then
                                if(abs(global_stiff1(iseq(totveq(itotv)))).le.1.e-15)  &
                                    global_stiff1(iseq(totveq(itotv)))=1.e30
                            endif
                        enddo      !20220618
                    endif !20240312 YL


                    if(type_load=='LOAD2'.or.(kstat==2.and.iiter.le.2).or.(type_load/='LOAD2'.and.kstat/=2.and.iiter==1).or.  &
                        (ngaps/=0.and.istatec==0)) then	  !! for temperature 20130510
                        call gpvar2_initial
                        if (ninit/=0.and.(kinit==2.and.iincs==1.and.istep==1)) then  !20201203
                            call eload_initialize
                            call eload_initial_stress
                            call force_release
                            where(totveq==0)
                                torel=0.0
                            endwhere
                        endif

                        call eload_initialize
                        call residu_f
                        call eload_field
                        if (uwcpl/=0) call eload_couple
                        if(nbspring>0) &
                            call eload_back_spring
                        if(ground_inf/=0)call semi_inf_load

                        if (stabpw==1)call stabload
                        call force_internal

                    endif    !! end for if(type_load=='LOAD2'...
                    if (mdiv/=1) then
                        if (idiv==1.and.iiter==1.and.allocated(torel))toform=toform+torel
                    else
                        if (iiter==1.and.allocated(torel))tofor=tofor+torel
                    endif

                    if(ngaps/=0.and.iblks>=iblks_bt.and.mdiv==1)call ctfor_to_tofor(tofor0,tofor)
                    if(ngaps/=0.and.iblks>=iblks_bt.and.mdiv/=1)call ctfor_to_tofor(tofor0,toform)

                    if(nrcsteel/=0.and.mdiv==1)call csfor_to_tofor(tofor0,tofor)
                    if(nrcsteel/=0.and.mdiv/=1)call csfor_to_tofor(tofor0,toform)

                    if(neuman==1.and.(istep/=1.or.iiter/=1).and.kresl/=0)  goto 2  !ctt2005 , change position!
                    if(type_nl==8.and.(iiter>1.or.(kstat==2.and.iiter>2))) goto 2  !MNR
                    if ((type_solver=='PROFILE'.or.type_solver=='PARDISO').and.kresl/=0)then
                        operation='FACTORIZE'
                        call solve

                        !20230216 形成反演边界温度需要的C矩阵
                        if(nbackdT==2.and.istep==1.and.iiter==1) then
                            if(allocated(cmatrix_dtv))deallocate(cmatrix_dtv)
                            if(allocated(inv_cmatrix_dtv2))deallocate(inv_cmatrix_dtv2)
                            allocate(cmatrix_dtv(Npoints_pbx,nfixsets),inv_cmatrix_dtv2(nfixsets,nfixsets))
                            call cmatrix_dtv_formation(cmatrix_dtv,inv_cmatrix_dtv2)
                        endif
                        !20230216



                    end if
2                   continue

                    logx=ngaps/=0.and.(iiter==1.and.istep==inc_step).and.iincs==(lincs+1) !20200331
                    if (logx)then !ctt2005
                        if (restart_ctt==0)then !restart_ctt
                            kdimn=ndimn
                            if(block_stab==1)kdimn=3*(ndimn-1) !2015/11/17
                            allocate(rot(kdimn,kdimn))
                            rot=0.
                            do igapb=1,ngapb
                                if(block_appear_process(igapb,iblks)==0)cycle    !20200331
                                npgblock=gapb(igapb)%npgblock
                                gapb(igapb)%cmatrix=0.

                                if(gapb(igapb)%eblock==0) cycle  !2017/11/19
                                do ipoin=1,npgblock
                                    igaps=gapb(igapb)%nodegblock_igaps(ipoin)
                                    ipair=gapb(igapb)%nodegblock_ipairs(ipoin)
                                    ij=gapb(igapb)%nodegblock_onetwo(ipoin)


                                    coef=1.
                                    if(ij==2)coef=-1.
                                    rot(1:ndimn,1:ndimn)=gaps(igaps)%rot(:,:,ipair)
                                    if(kdimn>ndimn)then
                                        if(ndimn==2)rot(3,3)=1.
                                        if(ndimn==3)rot(4:6,4:6)= rot(1:ndimn,1:ndimn)
                                    endif

                                    allocate(unitl(kdimn),unitg(kdimn))
                                    do idimn=1,kdimn

                                        itotvbt=(ipoin-1)*kdimn+idimn
                                        unitl=0.
                                        unitl(idimn)=1.*coef


                                        unitg=transpose(rot).x.unitl
                                        rvector=0.
                                        call  unit_force_trans(igapb,kdimn,ij,unitg,igaps,ipair,rvector)

                                        operation='SOLVE'
                                        call solve


                                        do kpoin=1,npgblock
                                            jgaps=gapb(igapb)%nodegblock_igaps(kpoin)
                                            jpair=gapb(igapb)%nodegblock_ipairs(kpoin)
                                            ij0=gapb(igapb)%nodegblock_onetwo(kpoin)         !2017/04/03
                                            call result_node_to_center(kdimn,ij0,jgaps,jpair,result,unitg)

                                            do jdimn=1,kdimn
                                                jtotvbt=(kpoin-1)*kdimn+jdimn
                                                gapb(igapb)%cmatrix(jtotvbt,itotvbt)=unitg(jdimn)
                                            end do
                                        end do  !kpoin
                                    end do  !idimn
                                    deallocate(unitl,unitg)
                                end do  !ipoin

                                do ipoin=1,npgblock
                                    igaps=gapb(igapb)%nodegblock_igaps(ipoin)
                                    ipair=gapb(igapb)%nodegblock_ipairs(ipoin)
                                    ij=gapb(igapb)%nodegblock_onetwo(ipoin)
                                    coef=1.
                                    if(ij==2)coef=-1.
                                    rot=0.
                                    rot(1:ndimn,1:ndimn)=gaps(igaps)%rot(:,:,ipair)

                                    if(kdimn>ndimn)then
                                        if(ndimn==2)rot(3,3)=1.
                                        !if(ndimn==3)rot(1:3,4:6)= rot(1:ndimn,1:ndimn)
                                        if(ndimn==3)rot(4:6,4:6)= rot(1:ndimn,1:ndimn)
                                        !if(ndimn==3)rot(4:6,1:3)= rot(1:ndimn,1:ndimn)
                                    endif

                                    rot=coef*rot

                                    allocate(cmatrixl(kdimn,kdimn))
                                    do jpoin=1,npgblock
                                        jtotv=(jpoin-1)*kdimn
                                        itotv=(ipoin-1)*kdimn
                                        cmatrixl=gapb(igapb)%cmatrix(itotv+1:itotv+kdimn,jtotv+1:jtotv+kdimn)
                                        gapb(igapb)%cmatrix(itotv+1:itotv+kdimn,jtotv+1:jtotv+kdimn)=   &
                                            rot.x.cmatrixl
                                    end do
                                    deallocate(cmatrixl)
                                end do
                            end do  !igapb
                            call forAdirect !fzx !形成A矩阵

                            do igapb=1,ngapb
                                !write(7,*)'igapb=',igapb,'ntotv_bt=',gapb(igapb)%ntotv_bt,'camatrix='
                                do itotvbt=1,gapb(igapb)%ntotv_bt
                                    !write(7,*)gapb(igapb)%cmatrix(itotvbt,:)
                                    do jtotvbt=1,gapb(igapb)%ntotv_bt
                                        write(recttunit)gapb(igapb)%cmatrix(itotvbt,jtotvbt)
                                    enddo
                                enddo
                            enddo  !igapb

                            deallocate(rot)

                        elseif(restart_ctt==1)then !restart_ctt
                            call forAdirect !fzx !形成A矩阵
                            rewind(recttunit)
                            do igapb=1,ngapb
                                npgblock=gapb(igapb)%npgblock
                                do itotvbt=1,gapb(igapb)%ntotv_bt
                                    do jtotvbt=1,gapb(igapb)%ntotv_bt
                                        read(recttunit)gapb(igapb)%cmatrix(itotvbt,jtotvbt)
                                    enddo
                                enddo
                            enddo
                        else !restart_ctt
                            write(*,*)'no such restart_ctt!!'
                            stop
                        endif !restart_ctt
                    endif  !!ctt2005

                    if(kresl/=0) & !20220311
                        call cmatrix_c_formation !20220311


                    rvector=0.0
                    if (type_solver/='JPCG') then
                        do itotv=1,ntotv
                            if (totveq(itotv)/=0)then
                                if (mdiv/=1) then
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                        toform(itotv)-stfor(itotv)
                                else
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                        tofor(itotv)-stfor(itotv)
                                endif
                                if(nonsym==0)then !20240312 YL
                                    if(abs(global_stiff1(iseq(totveq(itotv)))).le.1.e-15) &  !20220618
                                        global_stiff1(iseq(totveq(itotv)))=1.e20
                                endif !20240312 YL
                            endif
                        end do
                        if (nflow/=0)then
                            call flow_charge
                            do ii=1,nfreeflownode
                                ipoin=listfreeflownode(ii)
                                itotv=nodfn(lmdofn(8),ipoin)
                                if (totveq(itotv)/=0.and.(result_zero(itotv)>.1  &
                                    .or.(result_zero(itotv)>0..and.flowrate(ipoin)>0.))) then
                                    global_stiff1(iseq(totveq(itotv)))=1.e20
                                    rvector(totveq(itotv))=0.
                                    result_zero(itotv)=0.
                                endif
                            end do
                        endif
                        !!int2000
                        do itotv=1,ntotv
                            nintf=trans(itotv)%nintf
                            if (nintf/=0) then
                                iieq=totveq(itotv)
                                if (iieq/=0)rvector(iieq)=0.
                                do iintf=1,nintf
                                    iieq=totveq(trans(itotv)%listf(iintf))
                                    if (iieq/=0) &
                                        rvector(iieq)=rvector(iieq)+  &
                                        (tofor(itotv)-stfor(itotv))*trans(itotv)%rintf(iintf)
                                end do
                            endif
                        end do
                        !!int2000
                    else
                        if (mdiv/=1)rvector=toform-stfor
                        if (mdiv==1)rvector=tofor -stfor
                    endif

                    if (type_load=='ARCLENGTH'.and.kresl/=0) then
                        allocate(rvectorm(neq))
                        rvectorm=rvector
                        rvector=0.0
                        if (type_solver/='JPCG') then
                            do itotv=1,ntotv
                                if(totveq(itotv)/=0) &
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+tofor_arclength(itotv)
                            end do
                        else
                            rvector=tofor_arclength
                        endif
                        operation='SOLVE'
                        call solve
                        delta_arclength=result
                        rvector=rvectorm
                        deallocate(rvectorm)
                    endif

                    !write(7,*)'idiv=',idiv,'iiter=',iiter
                    !write(7,*)'rvector=',rvector

                    if (type_nl==8)then
                        if(kstat/=2)call bfgsr(iiter)
                        if(kstat==2)call bfgsr(iiter-1)
                    else
                        operation='SOLVE'
                        call solve
                    endif

                    !write(7,*)'result=',result

                    if(ngaps/=0.and.iblks>=abs(iblks_bt))call solve_ctt  !!ctt2005
                    if(nrcsteel/=0) call solve_bond_force_of_cs !20210328


                    if (type_load=='ARCLENGTH')then
                        call find_dfact_of_arclength (irst)
                        if (irst==1) then
                            time_begin=tcurves(arc_curve)%time_begin
                            detal=tcurves(arc_curve)%detal
                            if (abs(ttime-time_begin-ditime).le.1.e-8)detal=tcurves(arc_curve)%fact_inc
                            tcurves(arc_curve)%detal=detal*.5
                            if (abs(ttime-time_begin-ditime).le.1.e-8)tcurves(arc_curve)%fact_inc=detal*.5
                            goto 222
                        endif
                    endif

                    if(neuman==1.and.(istep/=1.or.iiter/=1).and.kresl/=0)call neuman_expan

                    call TIME(char_time)
                    print *, 'time: ', char_time
                    write(chkunit,*)'time: ', char_time



                    call varupdate

                    !20230216 反演边界温度
                    if(nbackdT==2.and.iiter==1) then
                        allocate(observstar(Npoints_pbx),dtv(nfixsets),dtvi(nfixsets))
                        observstar=0.

                        do i=1,mvalue
                            print *,i,'Value_observ(i)%istep=',Value_observ(i)%istep
                            if(Value_observ(i)%iblks/=iblks)cycle
                            if(Value_observ(i)%iincs/=iincs)cycle
                            if(Value_observ(i)%istep/=istep)cycle
                            ivalue_point=Value_observ(i)%ivalue_point
                            observstar(ivalue_point)=Value_observ(i)%value_measure
                            nintf=para_points(ivalue_point)%nintf
                            listf=>para_points(ivalue_point)%listf
                            rintf=>para_points(ivalue_point)%rintf
                            Tpredict=dot_product(result_zero(listf),rintf)
                            observstar(ivalue_point)=observstar(ivalue_point)-Tpredict
                            nullify(listf,rintf)
                        end do
                        print *,'observstar=',observstar,'Tpre=',Tpredict
                        dTv=transpose(cmatrix_dtv).x.observstar
                        dtvi=inv_cmatrix_dtv2.x.dtv


                        do i=1,ndofix    !20231130
                            mfixset=prescrib(i)%mfixset
                            do j=1,mfixset !20231130
                                ifixset=prescrib(i)%mlist(j)
                                idofn=prescrib(i)%ldofix
                                result_zero(idofn)= result_zero(idofn)+dtvi(ifixset)*prescrib(i)%rintf(j)
                            end do !20231130
                        end do !20231130

                        !   do i=1,ndofix
                        !ifixset=prescrib(i)%ifixset
                        !idofn=prescrib(i)%ldofix
                        !result_zero(idofn)= result_zero(idofn)+dtvi(ifixset)
                        !  end do

                        deallocate(observstar,dtv,dtvi)

                    endif
                    !20230216



                    call relative_dis_watertight !20231007 止水 !20240305
                    call eload_initialize
                    if(ikindks/=0) call strain_for_steel_bar !steel 2008
                    call residu_f

                    if(type_load/='LOAD2'.or.(type_load=='LOAD2'.and.idiv==2))then
                        call eload_field
                        if (uwcpl/=0) call eload_couple
                        if (stabpw==1)call stabload

                        call reaction_prescribed
                        call conver_load

                        if (nchek==0)call conver_nodal_value

                        if (nchek==0)exit
                    endif   !if(type_load/='LOAD2'...
                    if(type_load=='LOAD2'.and.idiv==1)goto 10
                end do   !! loop for iiter
10              continue

                if(type_load/='LOAD2'.or.(type_load=='LOAD2'.and.idiv==2))then   !20220607
                    call local_stress
                    call contact_state(1)
                    call crack_state !crack 2006
                    if(istatec==0) &
                        call state_and_stiff_2021
                    do igaps=1,ngaps
                        npairs=gaps(igaps)%npairs
                        do ipairs=1,npairs
                            if(gaps(igaps)%pair_process(ipairs)==0)cycle  !20200331
                            gaps(igaps)%dxyz0(:,ipairs)=gaps(igaps)%dxyz(:,ipairs)
                            gaps(igaps)%ctforce0(:,ipairs)=gaps(igaps)%ctforce(:,ipairs)
                            gaps(igaps)%state0(ipairs)=gaps(igaps)%state(ipairs)
                            gaps(igaps)%damage0(ipairs)=gaps(igaps)%damage(ipairs)
                            gaps(igaps)%kxyz0(:,:,ipairs)=gaps(igaps)%kxyz(:,:,ipairs)
                        end do
                    end do
                    if (nflow/=0)call flow_charge
                    call gpvarupdate
                endif   !20220607

            end do    !! for idiv

            if(modf_dis_blocks(iblks)==1)call construction_dis_modify

100         toforl=tofor
            if (istep/noutn*noutn==istep)then
                iwriten=iwriten+1
                call out_record
                call outputres !for output
            endif

            if(outinp<0)then  !稳定渗流场分析时向oip文件输出结点压力 20220623
                !write(outinpunit,'(a)')'ipoin    pore_pressure'
                !idofn=lmdofn(8)
                !do ipoin=1,npoin
                !    itotv=nodfn(idofn,ipoin)
                !    if(itotv==0)then
                !        write(outinpunit,'(i8,e16.6)')ipoin,0.0
                !    else
                !        pvalue=result_zero(itotv)
                !        if(pvalue<0)pvalue=0.0
                !        write(outinpunit,'(i8,e16.6)')ipoin,pvalue
                !    endif
                !end do
                !!!!
                allocate(midt(npoin))  !20220626
                icdofn=lmdofn(8)
                midt=0.
                do ipoin=1,npoin
                    itotv=nodfn(icdofn,ipoin)
                    if (itotv/=0) then
                        midt(ipoin)=result_zero(itotv)
                    endif
                end do
                write(outinpunit)midt
                deallocate(midt)
            endif


            if(nforce/=0.or.ngaps/=0)call force_interface
            if(kstab/=0.)call safety_factor


            if (istep/noutf*noutf==istep)then
                if(nforce/=0.or.ngaps/=0)call write_force_interface
                call out_full_write
                if (outplot(1:3)=='GID')   call OUT_GID_WRITE
                if (outplot(1:6)=='COSMOS')call OUT_COSMOS_WRITE
            endif

            if (istep/nresta*nresta==istep)call resta_read_write(-1)


            if(Bparameter>0)then !20230523
                !Value_observ(:)%value_computation=0.
                do ivalue=1,mvalue
                    !if(Value_observ(ivalue)%ic==0)cycle
                    iblks_i=Value_observ(ivalue)%iblks
                    iincs_i=Value_observ(ivalue)%iincs
                    istep_i=Value_observ(ivalue)%istep
                    idofn =lmdofn(Value_observ(ivalue)%idofn)
                    ivalue_point=Value_observ(ivalue)%ivalue_point
                    if(iblks_i==iblks.and.iincs_i==iincs.and.istep_i==istep)then
                        nintf=para_points(ivalue_point)%nintf
                        listf=>para_points(ivalue_point)%listf
                        rintf=>para_points(ivalue_point)%rintf
                        if(Bparameter==1)Value_observ(ivalue)%value_computation=dot_product(rintf,result_zero(nodfn(idofn,listf)))
                        if(Bparameter==2)Value_observ(ivalue)%value_computation=dot_product(rintf,deltafi(listf))
                        nullify(listf,rintf)
                    endif
                end do
            endif   !20230523

            if(Bparameter<0)then !20200812
                tbstep=tbstep+1
                do i=1,nback_point

                    bblks=freedom_for_back(4,i)  !20230523
                    if(bblks>iblks)cycle !20230523
                    inode=freedom_for_back(1,i)
                    idofn=freedom_for_back(2,i)
                    jnode=freedom_for_back(3,i)
                    print *,'inode=',inode,'idofn=',idofn,'jnode=',jnode
                    itotv=nodfn(lmdofn(idofn),inode)
                    if(jnode/=0)jtotv=nodfn(lmdofn(idofn),jnode)
                    if(Bparameter==-1)then
                        Value_vc(i,tbstep,istoch)=result_zero(itotv)
                        if(jnode/=0)Value_vc(i,tbstep,istoch)=Value_vc(i,tbstep,istoch)-result_zero(jtotv)
                    elseif(Bparameter==-2)then
                        Value_vc(i,tbstep,istoch)=deltafi(itotv)
                        if(jnode/=0)Value_vc(i,tbstep,istoch)=Value_vc(i,tbstep,istoch)-deltafi(jtotv)
                    end if
                end do
            endif  !20200812
            if((bparameter>=1.and.bparameter<=2).and.balgor>=1) call dudx

            if(upliftin<0)then  !考虑渗流场影响分析时向upf文件输出结点压力 20220623
                allocate(midt(npoin))  !20220626
                icdofn=lmdofn(8)
                midt=0.
                do ipoin=1,npoin
                    itotv=nodfn(icdofn,ipoin)
                    if (itotv/=0) then
                        midt(ipoin)=result_zero(itotv)
                    endif
                end do
                write(upliftunit)midt
                deallocate(midt)

            endif   !20230708


        end do     !! loop for istep

        if(Qstatic/=0) then !20221104
            deallocate(qstatic_force%appearg,qstatic_force%qfactor,qstatic_force%cor_coef) !20221104
            deallocate(qstatic_force)  !20221104
        endif !20221104
        if(cwater/=0.and.delgroup>0)deallocate(coef_water)
    end do     !! loop for iincs


    !if(winit==-1) call out_next_write
    if(winit==-1*iblks) call out_next_write !20231215YULI

    END SUBROUTINE STATIC_U_PW

    SUBROUTINE STATIC_U_PWm

    character(80)text
    integer(ink) itotv,ielem
    real   (irk) time
    integer(ink) iintf,nintf,iieq,idofn,ipoin   !!int2000

    print *,'static_U_Pwm'
    read(mainunit,*)text
    read(mainunit,*)nincs

    mdiv=2
    time=0.0
    do iincs=1,lincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta
        read(mainunit,*)toler_force,toler_var(1:mdofn)
    end do


    do iincs=lincs+1,nincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta
        read(mainunit,*)toler_force,toler_var(1:mdofn)

        do istep=1,nstep
            if (outintr>0.and.iblks>=outintr)trstep=trstep+1  !20200226
            time=time+ditime
            ttime=ttime+ditime        !! only for output

            call dfact_time_curve(ttime)

            call modf_var_prescribed
            call gravity
            call loadfl
            call force_external

            do idiv=1,mdiv

                deltafi=0.0
                delitfi=0.0

                if (mdiv/=1) &
                    toform=toforl+(tofor-toforl)*idiv/mdiv

                do iiter=1,miter

                    print *,'iincs=',iincs,'istep=',istep,'iiter=',iiter


                    call algort
                    if (iiter==1)then
                        deltafi=0.0
                        delitfi=0.0
                        call predict
                        do ielem=1,nelem   !!simo_rifai
                            if (associated(element(ielem)%alfa))element(ielem)%alfa=0.
                        end do  !!simo_rifai
                    endif
                    if (kresl/=0.and.uwcpl/=0)call stiff_u
                    if (ksmat/=0.and.uwcpl==1)call mcmatrx('W')
                    if (kmass/=0)             call mcmatrx('U')
                    if (khmat/=0.and.uwcpl/=1)call hmatrx('W')
                    if (kqmat/=0.and.uwcpl/=0)call upwcouple
                    if (stabpw==1.and.khmat/=0)call stabpatch    !!stablize



                    if (kresl/=0.or.ksmat/=0.or.khmat/=0.or.kqmat/=0) then
                        if (type_solver/='JPCG')global_stiff1=0.0
                        if (nonsym/=0.and.type_solver=='PROFILE')global_stiff2=0.0
                        if (type_solver=='JPCG')then
                            do ielem=1,nelem
                                element(ielem)%estif=0.0
                            end do
                        endif
                        call estif_assemble
                        if (uwcpl==1)  call couple_assemble
                        if (stabpw==1) call stabpw_assemble           !! stablize
                    endif


                    if (idiv==1.and.iiter==1.and.istep==inc_step) then    !! iiter==1 and istep==inc_step
                        call eload_initialize
                        if ((ninit/=0.and.((iblks==1.and.kinit==1)).or.(kinit==2.and.iincs==1))) then

                            call eload_initial_stress
                            call gpvar2_initial
                            if (kinit==2.and.iincs==1)then
                                call force_release
                                where(totveq==0)
                                    torel=0.0
                                endwhere
                            endif

                        endif

                        if (uwcpl/=0.and.(iblks/=1.or.(iblks==1.and.iincs/=1)))call residu_f
                        call eload_field
                        if (uwcpl/=0) call eload_couple
                        if (stabpw==1)call stabload
                        call force_internal

                    endif    !! end for iiter==1 and istep==inc_step
                    if (mdiv/=1) then
                        if (idiv==1.and.iiter==1.and.allocated(torel))toform=toform+torel
                    else
                        if (iiter==1.and.allocated(torel))tofor=tofor+torel
                    endif

                    rvector=0.0
                    if (type_solver/='JPCG') then
                        do itotv=1,ntotv
                            if (totveq(itotv)/=0)then
                                if (mdiv/=1) then
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                        toform(itotv)-stfor(itotv)
                                else
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                        tofor(itotv)-stfor(itotv)
                                endif
                            endif
                        end do

                        !!int2000
                        do itotv=1,ntotv
                            nintf=trans(itotv)%nintf
                            if (nintf/=0) then
                                iieq=totveq(itotv)
                                if (iieq/=0)rvector(iieq)=0.
                                do iintf=1,nintf
                                    iieq=totveq(trans(itotv)%listf(iintf))
                                    if (iieq/=0) &
                                        rvector(iieq)=rvector(iieq)+  &
                                        (tofor(itotv)-stfor(itotv))*trans(itotv)%rintf(iintf)
                                end do
                            endif
                        end do
                        !!int2000
                    else
                        if (mdiv/=1)rvector=toform-stfor
                        if (mdiv==1)rvector=tofor -stfor
                    endif

                    if (type_nl==8.and.(iiter>1.or.(kstat==2.and.iiter>2))) goto 2  !MNR
                    if (type_solver=='PROFILE'.and.               &
                        (kresl/=0.or.ksmat/=0.or.khmat/=0.or.kqmat/=0)) then
                        operation='FACTORIZE'
                        call solve
                    end if
2                   if(type_nl==8) then
                        if (kstat/=2)call bfgsr(iiter)
                        if (kstat==2)call bfgsr(iiter-1)
                    else
                        operation='SOLVE'
                        call solve
                    endif

                    if (idiv==1) then
                        do idofn=1,cdofn
                            do ipoin=1,npoin
                                itotv=nodfn(idofn,ipoin)
                                if (itotv/=0.and.iffix(itotv)==0)  then
                                    delitfi(itotv)=result(itotv)
                                    deltafi(itotv)=deltafi(itotv)+delitfi(itotv)
                                endif
                            enddo     !! for ipoin
                        end do     !! for idofn
                        call residu_f
                        goto 10
                    endif

                    call varupdate
                    call eload_initialize
                    if (uwcpl/=0)call residu_f
                    call eload_field
                    if (uwcpl/=0)call eload_couple
                    if (stabpw==1)call stabload

                    call reaction_prescribed
                    call conver_load

                    if (nchek==0)call conver_nodal_value
                    if (nchek==0)exit
                    call gpvarupdate
                end do   !! loop for iiter
10              continue
            end do    !! for idiv
100         toforl=tofor
            if (istep/noutn*noutn==istep)then
                iwriten=iwriten+1
                call out_record
            endif

            if (istep/noutf*noutf==istep)then
                call out_full_write
                if (outplot(1:3)=='GID')   call OUT_GID_WRITE
                if (outplot(1:6)=='COSMOS')call OUT_COSMOS_WRITE
            endif

            if (istep/nresta*nresta==istep)call resta_read_write(-1)
        end do     !! loop for istep

    end do     !! loop for iincs


    !if(winit==-1) call out_next_write
    if(winit==-1*iblks) call out_next_write !20231215YULI

    END SUBROUTINE STATIC_U_PWm
