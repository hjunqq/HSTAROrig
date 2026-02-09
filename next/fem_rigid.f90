    SUBROUTINE STATIC_rigid_reli  !20211016

    logical logx
    character(80)text
    integer(ink) itotv,ielem,irst,trstep0,ipoin,idofn,ij,idofix,ldofix,idelgroup,i0,ipairs
    integer(ink) ncmat,ncpld,nstoch
    real   (irk) xtime,time_begin,detal,ttime0,coef
    real   (irk),allocatable::rvectorm(:),value(:)
    integer(ink) iintf,nintf,iieq,ivcoh,ivfri   !!int2000

    integer(ink) igapb,npgblock,jpoin,igaps,ipair,idimn,itotvbt,jdimn, &  !!ctt2005
        jtotv,kpoin,lpoin,jtotvbt,npairs,cwater   !!ctt2005
    real   (irk),allocatable::rot(:,:),tofor0(:)  !!ctt2005
    real   (irk),allocatable::unitl(:),unitg(:),cmatrixl(:,:) !!ctt2005

    integer(ink) iter,mkiter,i,j0,nbeta,ibeta,ic,Nv
    real(irk)    beta,er,gy
    integer(ink),allocatable::ja(:)
    real(irk),allocatable::xa(:),ya(:),sd(:),ed(:),ep(:),sp(:),ga(:),rc(:),ee(:),ss(:), &
        cov(:,:)

    read(stocunit,*)text
    read(stocunit,*)nbeta,nv,mkiter !number of radom,maximum iteration number

    allocate(betas(nbeta))
    do ibeta=1,nbeta
        allocate(betas(ibeta)%ga(nv))
        betas(ibeta)%ga=0.
    end do
    allocate(xa(nv),ya(nv),ed(nv),sd(nv),ep(nv),sp(nv),ga(nv),rc(nv),ja(nv),ee(nv),ss(nv))
    allocate(cov(nv,nv))

    read(stocunit,*)text
    read(stocunit,*)ja(:)  !distribution type:1,normal,2,log normal,3,extreme
    read(stocunit,*)ee(:)  !average value
    read(stocunit,*)ss(:)  !variance
    do i=1,nv
        read(stocunit,*)cov(i,:)  !correlation matrix
    end do

    !if(meshc==1.or.rmesh/=0)rewind(mainunit)

    read(mainunit,*)text
    read(mainunit,*)nincs

    print *,' in static_rigid_1**'

    if(ngaps/=0)allocate(tofor0(ntotv)) !!ctt2005

    do iincs=1,lincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta,cwater
        read(mainunit,*)toler_force,toler_var(1:mdofn)
    end do

    xtime=0.0
    do iincs=lincs+1,nincs
        print *,'iincs=',iincs

        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta,cwater
        read(mainunit,*)toler_force,toler_var(1:mdofn)
        if(cwater/=0.and.delgroup>0)then
            allocate(coef_water(delgroup,nstep))
            do idelgroup=1,delgroup
                read(mainunit,*)i0,coef_water(idelgroup,:)
            end do
        end if

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

222             call force_external
                !if(iblks/=1)mdiv=1   !5
                if(type_load=='LOAD2')mdiv=2  !!806
                do idiv=1,mdiv
                    !! temperature
                    if(type_load=='LOAD2'.and.idiv==2) goto 71
                    call load_of_creep_and_temperature
                    call creep_strain_of_rock_fill    !20130510
71                  if(mdiv/=1)toform=toforl+(tofor-toforl)*idiv/mdiv
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
                            delitfi=0.0
                            call predict
                        endif

                        if(ngaps/=0.and.iblks>=iblks_bt.and.iiter==1.and.mdiv==1)call ctfor_to_tofor(tofor0,tofor)  !!ctt2005
                        if(ngaps/=0.and.iblks>=iblks_bt.and.iiter==1.and.mdiv/=1)call ctfor_to_tofor(tofor0,toform)  !!ctt2005

                        logx=ngaps/=0.and.(iiter==1.and.istep==inc_step).and.iblks==iblks_bt
                        if (logx)then !ctt2005
                            if (restart_ctt==0)then !restart_ctt
                                do igapb=1,ngapb
                                    gapb(igapb)%cmatrix=0.
                                end do

                                call forAdirect !fzx !形成A矩阵
                                do igapb=1,ngapb
                                    do itotvbt=1,gapb(igapb)%ntotv_bt
                                        do jtotvbt=1,gapb(igapb)%ntotv_bt
                                            write(recttunit)gapb(igapb)%cmatrix(itotvbt,jtotvbt)
                                        enddo
                                    enddo
                                enddo  !igapb

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

                        if(ngaps/=0.and.iblks>=abs(iblks_bt))call solve_ctt_rigid  !!ctt2005


                        call TIME(char_time)
                        print *, 'time: ', char_time
                        write(chkunit,*)'time: ', char_time


                        call varupdate
                        call conver_nodal_value
                        if(nchek==0)exit !tcl
10                      continue
                        print *,'miter=',miter,'iiter=',iiter
                    end do   !! loop for iiter

                    if(istatec==0) &
                        call state_and_stiff_rigid_2021
                    do igaps=1,ngaps
                        npairs=gaps(igaps)%npairs
                        do ipairs=1,npairs
                            if(gaps(igaps)%pair_process(ipairs)==0)cycle  !20200331
                            gaps(igaps)%dxyz0(:,ipairs)=gaps(igaps)%dxyz(:,ipairs)
                        end do
                    end do

                end do    !! for idiv


100             toforl=tofor

                if (istep/noutn*noutn==istep)then
                    iwriten=iwriten+1
                    call out_record
                    call outputres !for output
                endif


                if(nforce/=0.or.ngaps/=0)call force_interface  !2017/11/19
                if(kstab/=0.)call safety_factor  !2017/11/19

                if (istep/noutf*noutf==istep)then
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
        if(cwater/=0.and.delgroup>0)deallocate(coef_water)
    end do !!iincs
    if(sysrelis>0) &
        call system_reliability(nbeta,nv)

    !if(winit==-1) call out_next_write
    if(winit==-1*iblks) call out_next_write !20231215YULI
    if(ngaps/=0)deallocate(tofor0)  !!ctt2005

    END SUBROUTINE STATIC_rigid_reli !20211016

    SUBROUTINE STATIC_rigid_1

    logical logx
    character(80)text
    integer(ink) itotv,ielem,irst,trstep0,ipoin,idofn,ij,idofix,ldofix,idelgroup,i0,ipairs
    real   (irk) xtime,time_begin,detal,ttime0,coef
    real   (irk),allocatable::rvectorm(:),value(:)
    integer(ink) iintf,nintf,iieq   !!int2000

    integer(ink) igapb,npgblock,jpoin,igaps,ipair,idimn,itotvbt,jdimn, &  !!ctt2005
        jtotv,kpoin,lpoin,jtotvbt,npairs,cwater,Qstatic   !!ctt2005
    real   (irk),allocatable::rot(:,:),tofor0(:)  !!ctt2005
    real   (irk),allocatable::unitl(:),unitg(:),cmatrixl(:,:) !!ctt2005

    !if(meshc==1.or.rmesh/=0)rewind(mainunit)

    read(mainunit,*)text
    read(mainunit,*)nincs

    print *,' in static_rigid_1**'

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

222         call force_external
            !if(iblks/=1)mdiv=1   !5
            if(type_load=='LOAD2')mdiv=2  !!806
            do idiv=1,mdiv
                !! temperature
                if(type_load=='LOAD2'.and.idiv==2) goto 71
                call load_of_creep_and_temperature
                call creep_strain_of_rock_fill    !20130510
71              if(mdiv/=1)toform=toforl+(tofor-toforl)*idiv/mdiv
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
                        delitfi=0.0
                        call predict
                    endif

                    if(ngaps/=0.and.iblks>=iblks_bt.and.iiter==1.and.mdiv==1)call ctfor_to_tofor(tofor0,tofor)  !!ctt2005
                    if(ngaps/=0.and.iblks>=iblks_bt.and.iiter==1.and.mdiv/=1)call ctfor_to_tofor(tofor0,toform)  !!ctt2005

                    logx=ngaps/=0.and.(iiter==1.and.istep==inc_step).and.iblks==iblks_bt
                    if (logx)then !ctt2005
                        if (restart_ctt==0)then !restart_ctt
                            do igapb=1,ngapb
                                gapb(igapb)%cmatrix=0.
                            end do

                            call forAdirect !fzx !形成A矩阵
                            do igapb=1,ngapb
                                do itotvbt=1,gapb(igapb)%ntotv_bt
                                    do jtotvbt=1,gapb(igapb)%ntotv_bt
                                        write(recttunit)gapb(igapb)%cmatrix(itotvbt,jtotvbt)
                                    enddo
                                enddo
                            enddo  !igapb

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
90                  format(10e12.5)

                    if(ngaps/=0.and.iblks>=abs(iblks_bt))call solve_ctt_rigid  !!ctt2005


                    call TIME(char_time)
                    print *, 'time: ', char_time
                    write(chkunit,*)'time: ', char_time


                    call varupdate
                    call conver_nodal_value
                    if(nchek==0)exit !tcl
10                  continue
                    print *,'miter=',miter,'iiter=',iiter
                end do   !! loop for iiter

                if(istatec==0) &
                    call state_and_stiff_rigid_2021
                do igaps=1,ngaps
                    npairs=gaps(igaps)%npairs
                    do ipairs=1,npairs
                        if(gaps(igaps)%pair_process(ipairs)==0)cycle  !20200331
                        gaps(igaps)%dxyz0(:,ipairs)=gaps(igaps)%dxyz(:,ipairs)
                    end do
                end do

            end do    !! for idiv


100         toforl=tofor

            if (istep/noutn*noutn==istep)then
                iwriten=iwriten+1
                call out_record
                call outputres !for output
            endif


            if(nforce/=0.or.ngaps/=0)call force_interface  !2017/11/19
            if(kstab/=0.)call safety_factor  !2017/11/19

            if (istep/noutf*noutf==istep)then
                if(nforce/=0.or.ngaps/=0)call write_force_interface
                call out_full_write
                if(outplot(1:3)=='GID')   call OUT_GID_WRITE
                if(outplot(1:6)=='COSMOS')call OUT_COSMOS_WRITE
            endif

            if(istep/nresta*nresta==istep)call resta_read_write(-1)

        end do     !! loop for istep
        if(cwater/=0.and.delgroup>0)deallocate(coef_water)
    end do !!iincs
    !if(winit==-1) call out_next_write
    if(winit==-1*iblks) call out_next_write !20231215YULI
    if(ngaps/=0)deallocate(tofor0)  !!ctt2005

    END SUBROUTINE STATIC_rigid_1
