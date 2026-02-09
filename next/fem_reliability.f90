    SUBROUTINE STATIC_U_reli

    logical logx
    character(80)text,material,criteria
    integer(ink) itotv,ielem,irst,trstep0,ipoin,idofn,ij,ij0,idofix,ldofix,idelgroup,i0,ipairs
    integer(ink) ncmat,ncpld,nstoch
    real   (irk) xtime,time_begin,detal,ttime0,coef
    real   (irk),allocatable::rvectorm(:),value(:)
    integer(ink) iintf,nintf,iieq,ivcoh,ivfri   !!int2000

    integer(ink) igapb,npgblock,jpoin,igaps,ipair,idimn,itotvbt,jdimn, &  !!ctt2005
        jtotv,kpoin,lpoin,jtotvbt,npairs,cwater,jgaps,jpair,kdimn   !!ctt2005
    real   (irk),allocatable::rot(:,:),tofor0(:)  !!ctt2005
    real   (irk),allocatable::unitl(:),unitg(:),cmatrixl(:,:) !!ctt2005

    integer(ink) iter,mkiter,i,j0,nbeta,ibeta,ic,Nv
    real(irk)    beta,er,gy
    integer(ink),allocatable::ja(:)
    real(irk),allocatable::xa(:),ya(:),sd(:),ed(:),ep(:),sp(:),ga(:),rc(:),ee(:),ss(:), &
        cov(:,:)

    !
    !open(stocunit,file='stoc.dat')

    read(stocunit,*)text
    read(stocunit,*)nbeta,nv,mkiter !number of radom,maximum iteration number

    allocate(betas(nbeta))

    do ibeta=1,nbeta
        allocate(betas(ibeta)%ga(nv))
        betas(ibeta)%ga=0.
    end do

    !   read(stocunit,*)text
    !read(stocunit,*)group(1:ngroup)%ivcoh
    !   read(stocunit,*)group(1:ngroup)%ivfri


    allocate(xa(nv),ya(nv),ed(nv),sd(nv),ep(nv),sp(nv),ga(nv),rc(nv),ja(nv),ee(nv),ss(nv))
    allocate(cov(nv,nv))

    read(stocunit,*)text
    read(stocunit,*)ja(:)  !distribution type:1,normal,2,log normal,3,extreme
    read(stocunit,*)ee(:)  !average value
    read(stocunit,*)ss(:)  !variance
    do i=1,nv
        read(stocunit,*)cov(i,:)  !correlation matrix
    end do

    read(mainunit,*)text
    read(mainunit,*)nincs

    print *,' in static_U**'

    if(ngaps/=0)allocate(tofor0(ntotv)) !!ctt2005

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
        print *,'iincs=',iincs

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


        do ibeta=1,nbeta !!!2018/01/10
            read(stocunit,*)text
            do igaps=1,ngaps
                read(stocunit,*)gaps(igaps)%ivcoh,gaps(igaps)%ivfri
            end do

            iter=0
            SP=ss
            EP=ee
            XA=EP
            YA=0.
            er=1.e-3



333         iter=iter+1
            CALL DANGLI(nv,EE,SS,XA,JA,SD,ED,SP,EP)

            do igaps=1,ngaps
                if(nforce_gaps_appear(igaps)==2.or.nforce_gaps_appear(igaps)==0)cycle
                ivcoh=gaps(igaps)%ivcoh
                ivfri=gaps(igaps)%ivfri
                if(ivcoh==0.and.ivfri==0)cycle
                npairs=gaps(igaps)%npairs
                do ipairs=1,npairs
                    if(ivfri>0) &
                        gaps(igaps)%frict(ipairs)=xa(ivfri)
                    if(ivcoh>0) &
                        gaps(igaps)%cohes(ipairs)=xa(ivcoh)
                end do
            end do


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

                call gravity
                if(rmesh>0)call gravity1
                if(rmesh>1)call gravity2
                write(7,*)'cwater=',cwater,'delgroup=',delgroup
                if(cwater/=0.and.delgroup/=0)call step_water_pressure  !2013/3/18

                write(7,*)'force_external_reli**'
222             call force_external
                !if(iblks/=1)mdiv=1   !5
                if(type_load=='LOAD2')mdiv=2  !!806
                do idiv=1,mdiv
                    !! temperature
                    if(type_load=='LOAD2'.and.idiv==2) goto 71
                    call load_of_creep_and_temperature
                    call creep_strain_of_rock_fill    !20130510
71                  if(mdiv/=1)toform=toforl+(tofor-toforl)*idiv/mdiv
                    if(type_load=='DISCONTROL')preact0=prescrib(1)%rdofix
                    if(ngaps/=0.and.mdiv==1)tofor0=tofor !!ctt2005
                    if(ngaps/=0.and.mdiv/=1)tofor0=toform !!ctt2005
                    deltafi=0.0
                    do igapb=1,ngapb !fzx  tcl
                        if(gapb(igapb)%nrdof==0)cycle
                        gapb(igapb)%rdisp_deltafi=0.
                    enddo


                    do iiter=1,miter
                        iccontact=0 !zhao 05/07/30
                        print *,'iblks=',iblks,'idiv=',idiv,'iiter=',iiter

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
                                if(ground_inf/=0) call semi_inf_space_assemble

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

1                       if(type_load=='LOAD2'.or.(kstat==2.and.iiter.le.2).or.(type_load/='LOAD2'.and.kstat/=2.and.iiter==1).or.  &
                            (ngaps/=0.and.istatec==0)) then	  !! for temperature 20130510
                            if(type_load=='LOAD2'.and.idiv==2)then  !20130510
                                deltafi=0.0
                                delitfi=0.0
                            endif
                            call gpvar2_initial
                            if (ninit/=0.and.(kinit==2.and.iincs==1)) then
                                call eload_initialize
                                call eload_initial_stress
                                if (kinit==2.and.iincs==1)then
                                    call force_release
                                    where(totveq==0)
                                        torel=0.0
                                    endwhere
                                endif
                            endif

                            !write(7,*)'deltafi_before conver_load(100-120)'
                            !        do itotv=100,120
                            !        write(7,*)itotv,deltafi(itotv)
                            !        end do


                            call eload_initialize
                            call residu_f

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
                        write(*,*)'   ' !很奇怪，有这一行的话，就不出错，没有的话，SOLVE中allocate(resultm(ntotv))这一行出错！


                        if(ngaps/=0.and.iblks>=iblks_bt.and.iiter==1.and.mdiv==1)call ctfor_to_tofor(tofor0,tofor)  !!ctt2005
                        if(ngaps/=0.and.iblks>=iblks_bt.and.iiter==1.and.mdiv/=1)call ctfor_to_tofor(tofor0,toform)  !!ctt2005

                        !if(ngaps/=0)call ctfor_to_tofor(tofor0,tofor)  !!ctt2005

                        if(neuman==1.and.(istep/=1.or.iiter/=1).and.kresl/=0)  goto 2  !ctt2005 , change position!
                        if(type_nl==8.and.(iiter>1.or.(kstat==2.and.iiter>2))) goto 2  !MNR
                        if ((type_solver=='PROFILE'.or.type_solver=='PARDISO').and.kresl/=0)then
                            operation='FACTORIZE'
                            call solve
                        end if
2                       continue

                        logx=ngaps/=0.and.(iiter==1.and.istep==inc_step).and.iblks==iblks_bt
                        if (logx)then !ctt2005
                            if (restart_ctt==0)then !restart_ctt

                                kdimn=ndimn
                                if(block_stab==1)kdimn=3*(ndimn-1) !2015/11/17
                                allocate(rot(kdimn,kdimn))
                                rot=0.
                                do igapb=1,ngapb
                                    npgblock=gapb(igapb)%npgblock
                                    write(7,*)'igapb=',igapb,'npgblock=',npgblock
                                    gapb(igapb)%cmatrix=0.

                                    if(gapb(igapb)%eblock==0) cycle  !2017/11/19
                                    do ipoin=1,npgblock
                                        !jpoin=gapb(igapb)%nodegblock(ipoin)
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
                                            !
                                            ! if(igapb==2.and.idimn==2) then
                                            !write(7,*)'igapb**=',igapb,'jpoin=',gapb(igapb)%nodegblock(ipoin),'idimn=',idimn,'rvector='
                                            !      do jdimn=1,neq
                                            !          if(abs(rvector(jdimn))>1.e-6)write(7,*)jdimn,rvector(jdimn)
                                            !      end do
                                            ! endif

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


90                      format(10e12.5)
                        rvector=0.0
                        !write(7,*)'iiter=',iiter
                        !write(7,*)'rvector,tofor,stfor'
                        if (type_solver/='JPCG') then
                            do itotv=1,ntotv
                                if (totveq(itotv)/=0)then
                                    if (mdiv/=1)then
                                        rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                            toform(itotv)-stfor(itotv)

                                    else
                                        rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                            tofor(itotv)-stfor(itotv)
                                        !if(abs(rvector(totveq(itotv)))>1.e-8)write(7,*)itotv,tofor(itotv),stfor(itotv)

                                    endif
                                endif
                            end do

                            !!int2000
                            do itotv=1,ntotv
                                nintf=trans(itotv)%nintf
                                if (nintf/=0) then
                                    iieq=totveq(itotv)
                                    if(iieq/=0)rvector(iieq)=0.
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

                        if (type_nl==8)then
                            if(kstat/=2)call bfgsr(iiter)
                            if(kstat==2)call bfgsr(iiter-1)
                        else
                            operation='SOLVE'
                            call solve
                        endif

                        if(ngaps/=0.and.iblks>=abs(iblks_bt))call solve_ctt  !!ctt2005

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
                        call eload_initialize

                        if(ikindks/=0) call strain_for_steel_bar !steel 2008
                        !write(7,*)'bbxx residu_f'
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


                            if(nchek==0)exit !tcl
                        endif !ep2010
10                      continue
                        print *,'miter=',miter,'iiter=',iiter
                    end do   !! loop for iiter

                    if(type_load/='LOAD2') then
                        if(istatec==0) &
                            call state_and_stiff_2021
                        !write(7,*)'icttstif_static_u=',icttstif
                        do igaps=1,ngaps
                            npairs=gaps(igaps)%npairs
                            do ipairs=1,npairs
                                if(gaps(igaps)%pair_process(ipairs)==0)cycle  !20200331
                                gaps(igaps)%dxyz0(:,ipairs)=gaps(igaps)%dxyz(:,ipairs)
                                gaps(igaps)%ctforce0(:,ipairs)=gaps(igaps)%ctforce(:,ipairs)
                            end do
                        end do
                    endif

                    if(type_load/='LOAD2') &
                        call gpvarupdate

                    if(rmesh>0.and.nelem1>0)call gpvarupdate1
                    if(rmesh>1.and.nelem2>0)call gpvarupdate2
                end do    !! for idiv
                if(type_load=='LOAD2') &
                    call gpvarupdate
                if(modf_dis_blocks(iblks)==1)call construction_dis_modify


                if(type_load=='LOAD2') then
                    if(istatec==0) &
                        call state_and_stiff_2021
                    do igaps=1,ngaps
                        npairs=gaps(igaps)%npairs
                        do ipairs=1,npairs
                            if(gaps(igaps)%pair_process(ipairs)==0)cycle  !20200331
                            gaps(igaps)%state0(ipairs)=gaps(igaps)%state(ipairs)
                            gaps(igaps)%damage0(ipairs)=gaps(igaps)%damage(ipairs)
                            gaps(igaps)%ctforce0(:,ipairs)=gaps(igaps)%ctforce(:,ipairs)
                            gaps(igaps)%dxyz0(:,ipairs)=gaps(igaps)%dxyz(:,ipairs)
                        end do
                    end do
                endif


100             toforl=tofor

                if (istep/noutn*noutn==istep)then
                    iwriten=iwriten+1
                    call out_record
                    call outputres !for output
                endif
                !if (kstab==0.) then
                if(nforce/=0.or.ngaps/=0)call force_interface
                !else
                if(kstab/=0.)call safety_factor


                if (istep/noutf*noutf==istep)then
                    !if(kstab==0.and.nforce/=0)call write_force_interface
                    if(nforce/=0.or.ngaps/=0)call write_force_interface
                    call out_full_write
                    if(outplot(1:3)=='GID')   call OUT_GID_WRITE
                    if(outplot(1:6)=='COSMOS')call OUT_COSMOS_WRITE
                endif

                if(istep/nresta*nresta==istep)call resta_read_write(-1)
            end do     !! loop for istep

            call stab_rcandgy_reli(rc,xa,gy)
            call RI3(nv,GY,GA,RC,XA,YA,SD,ED,SP,EP,cov)  !translation between the distributions

            write(7,*)'iter=',iter,'gy=',gy,'xa=',xa
            if(abs(gy)>er.and.iter<mkiter) goto 333
            call betaindex(nv,rc,ep,sp,xa,cov,beta) ! computation reliability index
            write(7,*)'ibeta=',ibeta,'iter=',iter,'gy=',gy

            betas(ibeta)%beta=beta
            betas(ibeta)%ga(:)=-ga
            write(7,*)'reliability index=',beta
            write(7,*)'experiment points=',xa
            write(7,*)'alfa=',-ga


        end do !ibeta !!!2018/01/10
        !if(kstab==0.and.(nforce/=0.or.ngaps/=0))call write_force_interface
        if(Qstatic/=0) then !20221104
            deallocate(qstatic_force%appearg,qstatic_force%qfactor,qstatic_force%cor_coef) !20221104
            deallocate(qstatic_force)  !20221104
        endif !20221104

        if(cwater/=0.and.delgroup>0)deallocate(coef_water)
    end do !!iincs
    if(sysrelis>0) &
        call system_reliability(nbeta,nv)

    !if(winit==-1) call out_next_write
    if(winit==-1*iblks) call out_next_write !20231215YULI
    !if(ngaps/=0)deallocate(tofor0)  !!ctt2005

    END SUBROUTINE STATIC_U_reli

    SUBROUTINE betaindex(nv,rc,ep,sp,xa,cov,beta1)
    integer(ink) i,j,nv
    real(irk) sigmaz,miuz,beta1
    real(irk) rc(:),ep(:),sp(:),xa(:),cov(:,:)

    sigmaz=0.
    miuz=0.
    do i=1,nv
        miuz=miuz+(ep(i)-xa(i))*rc(i)
    end do
    do i=1,nv
        do j=1,nv
            sigmaz=sigmaz+rc(i)*sp(i)*rc(j)*sp(j)*cov(i,j)
        end do
    end do
    sigmaz=sqrt(sigmaz)
    beta1=miuz/sigmaz

    !	  write(2,*)'miuz=',miuz,'sigmaz=',sigmaz
    end SUBROUTINE

    SUBROUTINE RI3(nv,GY,GA,RC,XA,YA,SD,ED,SP,EP,cov)
    integer i,j,nv
    real(8) gg,y1,gy
    real(8) ga(:),rc(:),xa(:),ya(:),ed(:),sd(:),ep(:),sp(:),cov(:,:)
    real(8),allocatable::xc(:)
    allocate(xc(nv))
    xc=rc
    GG=0.0
    DO 160 I=1,NV
        rc(i)=0.
        do j=1,nv
            RC(I)=rc(i)+xc(j)*cov(i,j)*SP(j)
        end do
        GG=GG+RC(I)*RC(I)
160 CONTINUE
    GG=SQRT(GG)
    DO 164 I=1,NV
        IF(SP(I).EQ.0.OR.RC(I).EQ.0.) THEN
            GA(I)=0.
            YA(I)=0.
            GOTO 164
        ENDIF
        GA(I)=-RC(I)/GG
        YA(I)=SD(I)*YA(I)/SP(I)+(ED(I)-EP(I))/SP(I)
164 CONTINUE
    Y1=0.
    DO 165 I=1,NV
        Y1=Y1+YA(I)*GA(I)
165 CONTINUE
    Y1=Y1+GY/GG
    DO 170 I=1,NV
        YA(I)=GA(I)*Y1
170 CONTINUE
    DO 180 I=1,NV
        XA(I)=YA(I)*SP(I)+EP(I)
180 CONTINUE
    rc=xc

    deallocate(xc)
    END SUBROUTINE

    SUBROUTINE DANGLI(nv,EE,SS,XA,JA,SD,ED,SP,EP)
    integer(ink) i,j1,jj,nv
    real(irk) vv,sl,ar,ak,q,af,agf,aa,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10, &
        bb,yy,u1,uu,ub
    integer(ink) ja(:)
    real(irk) ee(:),ss(:),xa(:),sd(:),ed(:),ep(:),sp(:)
    DO 5 I=1,NV
        SD(I)=SP(I)
        ED(I)=EP(I)
        J1=JA(I)
        JJ=J1-2
        IF(JJ) 5,3,4
3       VV=(SS(I)/EE(I))**2    !
        SL=dLOG(1.+VV)
        IF(XA(I).GE.0.0) GO TO 6
        XA(I)=dABS(XA(I))
6       SP(I)=XA(I)*SQRT(SL)
        EP(I)=XA(I)*(1.+dLOG(EE(I))-dLOG(XA(I)*SQRT(1.+VV)))
        GOTO 5
4       AR=1.28255/SS(I)
        AK=EE(I)-0.5772/AR
        Q=EXP(-AR*(XA(I)-AK))
        AF=EXP(-Q)
        AGF=AR*Q*EXP(-Q)
        AA=1.-AF
        B0=1.570796288
        B1=3.706987906E-2
        B2=-8.364353589E-4
        B3=-2.250947176E-4
        B4=6.841218299E-6
        B5=5.824238515E-6
        B6=-1.04527497E-6
        B7=8.360937017E-8
        B8=-3.231081277E-9
        B9=3.657763036E-11
        B10=6.936233982E-13
        IF(AA-0.5) 230,240,250
230     BB=AA
        GO TO 255
250     BB=1.-AA
255     YY=-dLOG(4.*BB*(1.-BB))
        U1=YY*(B5+YY*(B6+YY*(B7+YY*(B8+YY*(B9+YY*B10)))))
        UU=YY*(B0+YY*(B1+YY*(B2+YY*(B3+YY*(B4+U1)))))
        UB=SQRT(UU)
        IF(AA.LT.0.5) GOTO 260
        UB=-UB
        GO TO 260
240     UB=0.0
260     SP(I)=EXP(-UB*UB/2.)/SQRT(2.*3.1416)/AGF
        EP(I)=XA(I)-SP(I)*UB
5   CONTINUE
    END SUBROUTINE

    subroutine system_reliability(nbeta,nv)
    integer(ink) nbeta,ntmod,itmod,jtmod,ibeta,i,j,nv
    real(irk) betax,pr,af,x0,lowpf,highpf,q1,q2,f0,relat_ave,beta_ave,aft
    integer(ink),allocatable:: rep(:)
    real(irk),allocatable::tga(:,:),tbeta(:),relat(:,:),ymult(:),zmult(:),qij(:,:)

    ntmod=nbeta
    allocate(tbeta(ntmod),tga(nv,ntmod),relat(ntmod,ntmod))
    relat=0.
    do ibeta=1,nbeta
        tga(:,ibeta)=betas(ibeta)%ga(:)
        tbeta(ibeta)=betas(ibeta)%beta
        write(7,*)'ibeta=',ibeta
        write(7,*)'ga=',tga(:,ibeta)
        write(7,*)'beta=',tbeta(ibeta)
    end do

    do itmod=1,ntmod
        do jtmod=1,ntmod
            relat(itmod,jtmod)=tga(:,itmod).d.tga(:,jtmod)
        end do
    end do
    write(7,*)' correlative matric'
    do itmod=1,ntmod
        write(7,10)relat(itmod,:)
    end do

    allocate(ymult(18),zmult(18))
    do i=1,18
        ymult(i)=1.
        do j=1,2*i+1,2
            ymult(i)=ymult(i)*real(j)
        end do
    end do
    zmult(1)=1.
    do i=2,18
        zmult(i)=1.
        do j=2*i-3,2*i-1,2
            zmult(i)=zmult(i)*real(j)
        end do
    end do
    !!!!for 串并联系统
    relat_ave=0.
    do i=1,ntmod
        do j=1,ntmod
            if(i/=j)relat_ave=relat_ave+relat(i,j)
        end do
    end do
    relat_ave=relat_ave/(ntmod*(ntmod-1))
    pr=1.0
    do i=1,ntmod
        !		  call af_normalx(tbeta(i),af,ymult,zmult)
        call af_normal(tbeta(i),af)
        pr=pr*af
    end do
    call af_beta(pr,beta_ave,ymult,zmult)


    call aft_gauss(beta_ave,relat_ave,aft,ntmod,ymult,zmult)
    !call aft_gauss(4.23_irk,0.57_irk,aft,10,ymult,zmult)
    print*,'串联系统失效概率=',aft
    write(7,*)'串联系统失效概率=',aft

    call af_beta(1-aft,betax,ymult,zmult)
    print*,'串联系统可靠指标=',betax
    write(7,*)'串联系统可靠指标=',betax

    betax=beta_ave*sqrt(ntmod/(1+relat_ave*(ntmod-1)))
    !betax=4.23*sqrt(10/(1+0.57*(10-1)))
    print*,'并联系统可靠指标=',betax
    write(7,*)'并联系统可靠指标=',betax
    !stop
    !!!!end for 串并联系统

    allocate(rep(ntmod))
    rep=0
    do i=1,ntmod
        if(rep(i)==0)then
            do j=i+1,ntmod
                if(rep(j)==0)then
                    if(relat(i,j)>0.7)then
                        rep(j)=1
                    endif
                endif
            end do
        endif
    end do
    pr=1.0
    do i=1,ntmod
        if(rep(i)==0)then
            !		  call af_normalx(tbeta(i),af,ymult,zmult)
            call af_normal(tbeta(i),af)
            pr=pr*af
        end if
    end do
    write(7,*)'PNET system reliability=',pr
    call af_beta(pr,betax,ymult,zmult)
    write(7,*)'PNET system reliability index=',betax
    !!!!following is for general narrow field method
    allocate(qij(ntmod,ntmod))
    do i=1, ntmod
        !			call af_normalx(tbeta(i),af,ymult,zmult)
        call af_normal(tbeta(i),af)
        qij(i,i)=1-af
    enddo
    do i=1, ntmod
        do j=i+1,ntmod
            x0=(tbeta(j)-relat(i,j)*tbeta(i))/(1-relat(i,j)**2)
            f0=1.
            if(x0<0.)then
                f0=-1.
                x0=x0*f0
            endif
            !			call af_normalx(x0,af,ymult,zmult)
            call af_normal(x0,af)
            write(7,*)'i=','j=',j,'x01=',x0,'af=',af
            if(f0<0.)af=1-af
            qij(i,j)=qij(i,i)*(1-af)
            x0=(tbeta(i)-relat(i,j)*tbeta(j))/(1-relat(i,j)**2)
            f0=1.
            if(x0<0.)then
                f0=-1.
                x0=x0*f0
            endif
            !			call af_normalx(x0,af,ymult,zmult)
            call af_normal(x0,af)
            write(7,*)'i=','j=',j,'x01=',x0,'af=',af
            if(f0<0.)af=1-af
            qij(j,i)=qij(j,j)*(1-af)
        end do
    end do
    write(7,*)'General narrow field method'
    write(7,*)'Qij****'
    do i=1,ntmod
        write(7,*)qij(i,:)
    end do

    lowpf=qij(1,1)
    do i=2,ntmod
        lowpf=lowpf+qij(i,i)
        do j=1,i-1
            lowpf=lowpf-qij(i,j)-qij(j,i)
        end do
    end do
    highpf=0.
    do i=1,ntmod
        highpf=highpf+qij(i,i)
    end do

    do i=2,ntmod
        q1=maxval(qij(1:(i-1),i))
        q2=maxval(qij(i,1:(i-1)))
        if(q2>q1)q1=q2
        highpf=highpf-q1
    end do
    write(7,*)'low  failure probability=',lowpf
    call af_beta(1-lowpf,betax,ymult,zmult)
    write(7,*)'high probability index=',betax
    write(7,*)'high failure probability=',highpf
    call af_beta(1-highpf,betax,ymult,zmult)
    write(7,*)'low probability index=',betax


    !!!!!!!end for


10  format(20e15.5)

    end subroutine system_reliability

    subroutine aft_gauss(beta_ave,relat_ave,aft,ntmod,ymult,zmult)
    integer(ink) i,ntmod
    real(irk) x0,y0,af ,beta_ave,relat_ave,aft,pi
    real(irk) w24(24),li24(24),ymult(:),zmult(:)
    pi=3.14159
    li24(1 )=-.06405689286260562609
    li24(2 )=-.19111886747361630916
    li24(3 )=-.31504267969616337439
    li24(4 )=-.43379350762604513849
    li24(5 )=-.54542147138883953566
    li24(6 )=-.64809365193697556925
    li24(7 )=-.74012419157855436424
    li24(8 )=-.82000198597390292195
    li24(9 )=-.88641552700440103421
    li24(10)=-.93827455200273275852
    li24(11)=-.97472855597130949820
    li24(12)=-.99518721999702136018
    li24(13)=.06405689286260562609
    li24(14)=.19111886747361630916
    li24(15)=.31504267969616337439
    li24(16)=.43379350762604513849
    li24(17)=.54542147138883953566
    li24(18)=.64809365193697556925
    li24(19)=.74012419157855436424
    li24(20)=.82000198597390292195
    li24(21)=.88641552700440103421
    li24(22)=.93827455200273275852
    li24(23)=.97472855597130949820
    li24(24)=.99518721999702136018
    w24(1 )=.12793819534675215697
    w24(2 )=.12583745634682829612
    w24(3 )=.12167047292780339120
    w24(4 )=.11550566805372560135
    w24(5 )=.10744427011596563478
    w24(6 )=.09761865210411388827
    w24(7 )=.08619016153195327592
    w24(8 )=.07334648141108030573
    w24(9 )=.05929858491543678075
    w24(10)=.04427743881741980617
    w24(11)=.02853138862893366318
    w24(12)=.01234122979998719955
    w24(13)=.12793819534675215697
    w24(14)=.12583745634682829612
    w24(15)=.12167047292780339120
    w24(16)=.11550566805372560135
    w24(17)=.10744427011596563478
    w24(18)=.09761865210411388827
    w24(19)=.08619016153195327592
    w24(20)=.07334648141108030573
    w24(21)=.05929858491543678075
    w24(22)=.04427743881741980617
    w24(23)=.02853138862893366318
    w24(24)=.01234122979998719955

    aft=1.
    do i=1,24
        x0=5*li24(i)
        y0=(beta_ave+sqrt(relat_ave)*li24(i)*5)/sqrt(1-relat_ave)
        !	  call af_normalx(y0,af,ymult,zmult)
        call af_normal(y0,af)
        af=5*af**ntmod*exp(-.5*x0**2)/sqrt(2.*pi)
        aft=aft-af*w24(i)
    end do


    end subroutine aft_gauss

    SUBROUTINE af_normalx(xy0,af,ymult,zmult)
    integer(ink) i,j
    real(irk) pi,af,xy0,y0,ymult(:),zmult(:)

    PI=3.141592653589793238
    if(xy0<=2.5)then
        y0=1.
        do j=1,18
            y0=y0+(xy0**(2.*j))/ymult(j)
        end do
        af=.5+xy0*exp(-.5*xy0**2)*y0/sqrt(2*pi)
    else
        y0=1.
        do j=1,18
            y0=y0+(-1)**j*(xy0**(-2.*j))*zmult(j)
        end do
        af=1.0-exp(-.5*xy0**2)*y0/xy0/sqrt(2*pi)
    endif

    END SUBROUTINE

    SUBROUTINE af_normal(xy0,af)
    integer i,j,k
    real(8) pi,af,xy0,y0

    pi=3.141592653589793238
    k=30
    if(xy0<=3.)then
        y0=(2.*k-1.)+k*(-1)**k*xy0**2/(2.*k+1)
        do j=k-1,1,-1
            y0=(2.*j-1.)+j*(-1)**j*xy0**2/y0
        end do
        af=.5+xy0*exp(-.5*xy0**2)/y0/sqrt(2*pi)
    else
        y0=xy0+k/xy0
        do j=k-1,1,-1
            y0=xy0+j/y0
        end do
        af=1-exp(-.5*xy0**2)/y0/sqrt(2*pi)
    endif

    END SUBROUTINE

    SUBROUTINE af_beta(af,beta,ymult,zmult)
    real(irk) pi,af,y,xy0,x1,af1,u0,beta,gx,af0,ymult(:),zmult(:)

    PI=3.141592653589793238

    af0=af
    if(af<0.5)af=1-af0

    y=-dlog(4*af*(1-af))
    x1=sqrt(y*(2.0611786-5.7262204/(y+11.640595)))
    if(x1<=4.5)then
        xy0=x1
    else
        xy0=x1-.000637*x1**2-.010437*x1+.059374
    endif
    !      call af_normalx(xy0,af1,ymult,zmult)
    call af_normal(xy0,af1)
    gx=af1-af
    u0=exp(-.5*xy0**2)/sqrt(2*pi)
    beta=xy0-gx/u0*(1-.5*xy0*gx/u0)
    if(af0<0.5)then
        beta=-beta
        af=af0
    endif

    END SUBROUTINE

    subroutine stab_rcandgy_reli(rc,xa,gy) !2018/01/10
    character(20)material,criteria
    integer(ink) ielem,igroup,iforce,jgroup,matno,ielgroup
    integer(ink) index,ngaus,igaus,order_int,nstre,npairs,igaps,ipairs,ivfri,ivcoh
    real   (irk),allocatable::ftang(:),fresi(:),ftang_gaps(:),fresi_gaps(:)
    real   (irk) ft,fn,uniax,frict,k_safety,elcod_local,yld,aera, dilan,sigma0,ftangt,fresit,t1,t2,tt,rc(:),xa(:),Gy

    rc=0.
    if(nforce==0.and.ngaps==0) return
    if(nforce==0) goto 10
    allocate(ftang(nforce),fresi(nforce))
    ftang=0.
    fresi=0.

    do iforce=1,nforce
        if(nforce_appear(iforce)/=1.and.nforce_appear(iforce)/=3)cycle !zhao 2010
        do jgroup=1,surface_force(iforce)%lgroup
            igroup=surface_force(iforce)%list(jgroup)
            index = group(igroup)%index
            nstre =group(igroup)%nstre
            matno =group(igroup)%matno
            order_int=elkn(index)%el_field(1)%order_intrules(1)
            ngaus = elkn(index)%ggaus(order_int)%ngaus
            elcod_local=group(igroup)%elcod_local
            uniax   =props(matno)%mechanical%solid%classicalEP%sigma0
            frict   =props(matno)%mechanical%solid%classicalEP%frict_angle
            ivcoh=group(igroup)%ivcoh
            ivfri=group(igroup)%ivfri
            if(ivcoh==0.or.ivfri==0)cycle
            DO ielgroup = 1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                do igaus=1,ngaus
                    yld=element(ielem)%field(1)%gpvar(nstre+2,igaus)
                    if (yld/=2.) then
                        aera=1./elcod_local*element(ielem)%egaus(1)%djacb(igaus)
                        fn=element(ielem)%field(1)%ntstress(1,igaus)
                        ft=element(ielem)%field(1)%ntstress(2,igaus)
                        fresi(iforce)=fresi(iforce)-fn*aera*xa(ivfri)+xa(ivcoh)*aera
                        ftang(iforce)=ftang(iforce)+ft*aera
                        rc(ivcoh)=rc(ivcoh)+aera
                        rc(ivfri)=rc(ivfri)-fn*aera

                    endif
                end do
            end do
        end do
        write(chkunit,*)'iforce=',iforce
        write(chkunit,*)'fresi=',fresi(iforce),'ftang=',ftang(iforce)
    end do !iforce

10  continue
    if(ngaps/=0)then
        allocate(fresi_gaps(ngaps),ftang_gaps(ngaps))
        fresi_gaps=0.;ftang_gaps=0.
        do igaps=1,ngaps
            if(nforce_gaps_appear(igaps)==2.or.nforce_gaps_appear(igaps)==0)cycle
            ivcoh=gaps(igaps)%ivcoh
            ivfri=gaps(igaps)%ivfri
            if(ivcoh==0.or.ivfri==0)cycle
            npairs=gaps(igaps)%npairs
            do ipairs=1,npairs

                fresi_gaps(igaps)=fresi_gaps(igaps)+gaps(igaps)%aera(ipairs)*xa(ivcoh)-gaps(igaps)%ctforce(ndimn,ipairs)*xa(ivfri)
                rc(ivcoh)=rc(ivcoh)+gaps(igaps)%aera(ipairs)
                rc(ivfri)=rc(ivfri)-gaps(igaps)%ctforce(ndimn,ipairs)
                t1=gaps(igaps)%ctforce(1,ipairs)
                if (ndimn==3)t2=gaps(igaps)%ctforce(2,ipairs)
                tt=abs(t1)
                if (ndimn==3)tt=sqrt(t1**2+t2**2)
                ftang_gaps(igaps)=ftang_gaps(igaps)+tt
            end do
        enddo
    end if

    fresit=0.;ftangt=0.
    if(nforce/=0)fresit=sum(fresi)
    if(ngaps/=0) fresit=fresit+sum(fresi_gaps)
    if(nforce/=0)ftangt=sum(ftang)
    if(ngaps/=0) ftangt=ftangt+sum(ftang_gaps)
    if(ngaps/=0)then
        write(chkunit,*)'fresit=',fresit,'ftangt=',ftangt
        write(chkunit,*)'fresi_gaps=',fresi_gaps
        write(chkunit,*)'ftang_gaps=',ftang_gaps
    endif

    Gy=fresit-ftangt

    write(7,*)'Gy=',Gy,'Rc=',Rc

    if(nforce/=0) &
        deallocate(ftang,fresi)
    if(ngaps/=0) &
        deallocate(ftang_gaps,fresi_gaps)

    end subroutine stab_rcandgy_reli  !2018/01/10
