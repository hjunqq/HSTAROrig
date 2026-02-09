    SUBROUTINE time_dependent

    character(80)text,type_curve
    integer(ink) itotv,ielem,itcurve,idofix,iextrf,nbounods,nnode,ibounod,ij,npairs,ipairs,ldofix,imcon,i0,ij0
    integer(ink) ilink,node1,node2,ipoin,idofn,kdimn,ivalue_point,bblks,nbasef,do_base_freq  !20230216
    real   (irk),allocatable::resultm(:),result0(:),tofor0(:)
    real   (irk),allocatable::disA(:),disB(:) !hxl2006 MIF
    real   (irk) time,coef,pvalue,Tpredict !20230216
    integer(ink) iintf,nintf,iieq,icdofn,ifixset   !!20230216
    integer(ink) igapb,npgblock,jpoin,igaps,ipair,jgaps,jpair,idimn,itotvbt,jdimn, &  !!ctt2005
        jtotv,kpoin,lpoin,jtotvbt,mfixset,j   !!20231130
    real   (irk),allocatable::rot(:,:)  !!ctt2005
    real   (irk),allocatable::unitl(:),unitg(:),cmatrixl(:,:),bounfreez(:),gaptofor(:),force_rigid(:)!!ctt2005
    integer(irk),allocatable::list_bound(:)
    real   (irk) s
    integer(ink) ilaymif
    integer(ink) i,iincs_i,iblks_i,istep_i,inode,jnode,ivalue  !20200819
    integer(ink),pointer::listf(:)  !20200819
    real   (irk),pointer::rintf(:)  !20200819
    real   (irk),allocatable::observstar(:),dtv(:),dtvi(:) !20230216

    integer(ink) cdtest,ntram,icrev,ICEND,Iseg,itram,icstop,nc1,nc2,nstre  !cdtest,2013/12/11
    real   (irk) p,q,eta,cvv,xl0,distf,dist,dist1,dist2,dist0,damage
    real   (irk),allocatable::xl1(:),xl2(:),midt(:) !20220626
    integer(ink),allocatable::ncyc(:),ic(:)

    real   (irk),allocatable::cmatrix_dtv(:,:),inv_cmatrix_dtv2(:,:)  !20230216

    integer(ink),pointer::ldofixb(:)  !hxl2006 MIF
    integer(ink),pointer::lnods(:)

    !if(outintw/=0)allocate(resultm(npoin),result0(ntotv))
    allocate(resultm(npoin),result0(ntotv))  !20200821
    if(allocated(earthquake_curve))deallocate(earthquake_curve)
    if(allocated(earthquake_curve_d))deallocate(earthquake_curve_d)
    if(allocated(earthquake_curve_v))deallocate(earthquake_curve_v)
    if(allocated(earthquake_curve_MIF))deallocate(earthquake_curve_MIF)
    if(allocated(fachv))deallocate(fachv)

    if (meshc==1.or.rmesh/=0)rewind(mainunit)
    if(Bparameter/=0.and.iblks==1)rewind(mainunit)  !20230902
    if(Bparameter/=0.and.iblks==1)rewind(upliftunit) !20230902


10  format(10i8)
    if(iblks==1)call cvoid(0)
    if(iblks==1.and.restart==0) then  !20220721
        if(any(group(:)%fieldid=='UW'))call cvoid(0)
        call inivdval
    endif

    allocate(earthquake_curve(ndimn),fachv(ndimn),earthquake_curve_d(ndimn),earthquake_curve_v(ndimn))!hxl2006 , VIE
    allocate(earthquake_curve_MIF(ndimn))
    earthquake_curve=0 ; earthquake_curve_MIF=0 ; earthquake_curve_d=0 ; earthquake_curve_v=0
    read(mainunit,*)text
    read(mainunit,*)nincs,cdtest,earthquake_curve(1:ndimn)
    read(mainunit,*)nstepjq,nstepjp !20231215YL

    if(type_ABC=='VIE')read(mainunit,*)hwdirec,hcoord  !20220105
    if(type_ABC=='VIE')read(mainunit,*)inpcord,earthquake_curve_d(1:ndimn),earthquake_curve_v(1:ndimn)
    if(type_ABC=='MIF')read(mainunit,*)inpcord,earthquake_curve_MIF(1:ndimn)

    mdiv=1
    if(cdtest>=1)then  !2013/12/11
        read(mainunit,*)text
        read(mainunit,*)ntram
        allocate(xl1(ntram),xl2(ntram),ncyc(ntram),ic(ntram))
        ic=0
        do itram=1,ntram
            read(mainunit,*)xl1(itram),xl2(itram)  !,ncyc(itram)
        end do
    endif

    !hxl2006 VIE
    allocate(freez(npoin))
    freez=0.0 ; nbounods=0
    if(type_ABC=='VIE')read(punit,*)nbounods
    write(7,*)'nbounods=',nbounods
    if (nbounods/=0) then
        allocate(list_bound(nbounods),bounfreez(nbounods))
        read(punit,*)list_bound(1:nbounods)
        read(punit,*)bounfreez(1:nbounods)

        do ibounod=1,nbounods
            ipoin=list_bound(ibounod)
            freez(ipoin)=bounfreez(ibounod)
        end do

        do ielem=1,nabssgroup
            lnods=>tabss(ielem)%lnods
            nnode=size(lnods)
            allocate(tabss(ielem)%cordzfree(nnode))
            do inode=1,nnode
                ipoin=lnods(inode)
                tabss(ielem)%cordzfree(inode)=freez(ipoin)
            end do
        end do
    end if
    deallocate(freez)
    !hxl2006 VIE


    !20231215YL
    if(nstepjq/=0)then  !20231008  file to sore the element average shear stress for judgement of liquifaction
        lquunit=72
        open(lquunit,file=probn(1:len1)//'.lqu',FORM='UNFORMATTED')
        nliqu=0
        disunit=73
        open(disunit,file=probn(1:len1)//'d.dis')
    endif  !
    if(nstepjp/=0)then
        pmtunit=74
        open(pmtunit,file=probn(1:len1)//'.pmt')
    endif
    omgunit=178 !20231008
    open(omgunit,file=probn(1:len1)//'.omg')
    !20231215YL

    if(cdtest>=1)then  !2013/12/11
        ICEND=0
        ICREV=0
        ISEG=0
        ITRAM=1
        ICSTOP=0
        fincre=1.
    endif



    if(ngaps/=0.or.nrcsteel/=0.or.nwcpipe/=0)allocate(tofor0(ntotv)) !!ctt2005

    time=0.0

    do iincs=1,lincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta,nmcon,nbasef  !20231215YL
        read(mainunit,*)toler_force,toler_var(1:mdofn)
    end do

    do iincs=lincs+1,nincs
        print *,'  time depend     iincs=',  iincs,'nincs=',nincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta,nmcon,nbasef !20231215YL
        read(mainunit,*)toler_force,toler_var(1:mdofn)
        if(nmcon/=0)allocate(lmcon(nmcon),rmcon(ndimn,nmcon))
        if (nmcon/=0) then
            read(mainunit,*)coef
            do imcon=1,nmcon
                read(mainunit,*)i0,lmcon(imcon),rmcon(:,imcon)
            end do
            rmcon=abs(rmcon)*coef
        endif


        !hxl2006 MIF

        !allocate(dissanru(nstep,ntotv),dissanzi(nstep,ntotv))   !sanshe  hxl
        allocate(disA(ntotv),disB(ntotv),disA_1(ntotv),disB_1(ntotv),disA_2(ntotv),disB_2(ntotv)) !hxl2006 MIF
        disA=0. ; disB=0. ; disA_1=0. ; disB_1=0. ; disA_2=0. ; disB_2=0.

        !dissanru=0.0
        !dissanzi=0.0                      !sanshe     hxl
        if (ntrans>0)then
            s=camif*ditime/dxmif
            k1(1)=(s-3)*s/2+1
            k1(2)=(2-s)*s
            k1(3)=(s-1)*s/2
            k2(1)=(2*s-3)*s+1
            k2(2)=4*(1-s)*s
            k2(3)=(2*s-1)*s
            k3(1)=4.5*(s-3.0)*s+1.
            k3(2)=3.0*(2.0-3.0*s)*s
            k3(3)=1.5*(3.0*s-1.0)*s
        endif
        !end hxl2006 MIF



        do istep=1,nstep

            print *, 'time: ', char_time
            !write(7,*)'gaps(1)%ctforce(:,1)=',gaps(1)%ctforce(:,1)
            disA_2=disA_1  !zhao
            disB_2=disB_1
            disA_1=disA
            disB_1=disB
            print *, 'iblks=',iblks,'iincs=',iincs,'istep=',istep
            if(outintr>0.and.iblks>=outintr)trstep=trstep+1  !20200226
            time=time+ditime
            ttime=ttime+ditime        !! only for output

            call dfact_time_curve(ttime)

            !!20231215YL
            do_base_freq=0  !20231008
            if(nbasef>0)then  !903
                if(istep==1.or.(istep/nbasef*nbasef==istep))then
                    if(gamamax==0)then !yuanli20230926
                        do_base_freq=1 !903
                        call base_frequency_analysis  !903
                    elseif(gamamax/=0.and.istep==1)then
                        do_base_freq=1 !903
                        call base_frequency_analysis  !903
                    endif
                endif
            elseif(nbasef<0)then  !903
                if(istep==1.or.(istep/abs(nbasef)*abs(nbasef)==istep))then
                    do_base_freq=1 !903
                    read(omgunit,*)text,base_freq
                    write(chkunit,*)'base_freq=',base_freq
                endif
            endif  !903
            !!20231215YL

            where(earthquake_curve==0)
                fachv=0.0
            elsewhere
                fachv=tcurves(earthquake_curve)%dfact
            endwhere
            write(7,*)'fachv=',fachv
            !levelset
            !***********************************************************************
            if (level_set_problem==2)then
                call levelsetmain
                ditime=delt
            endif
            !***********************************************************************
            if(outintr==0)then !20200220
                if(iincs==1.and.istep==inc_step)then
                    call placement_temperature(result0)  !20200220
                else
                    result0=result_zero
                endif
            endif
            allocate(inpru(ntotv),inpzi(ntotv))     !sanshe  hxl !hxl2006 MIF
            inpru=0.0 ; inpzi=0.0
            if(type_abc=='MIF')call modf_inpwav     !hxl2006 MIF
            !modf_inpwav：得到入射波场inpru及自由场inpzi

            deltafi=0.0
            call modf_var_prescribed
            if(submodel==1)call value_submodel_boundary  !20210321
            !modf_var_prescribed：在这个子程序里实现插值，由前几步人工边界区点的位移值得到当前步人工边界点的值，
            !以作为约束值，存在fixed中，或者说更新fixed
            !write(7,*)'result_zero(1:10)1=',result_zero(1:10)

            call heat_internal1

            do ielem=1,nelem   !!simo_rifai
                if(associated(element(ielem)%alfa))element(ielem)%alfa=0.
            end do  !!simo_rifai
            iccontact=0 !zhao 05/07/22

            do iiter=1,miter
                print *,'iiter=',iiter

                call algort
                if (istep==1.and.iiter==1)then
                    call local_stress
                    call contact_state(0)
                endif

                call porepr
                if(iiter==1)   call effect_stres_modul_for_steel_beam !20211125
                if(iiter==1)   call stiffness_for_bolt_spring  !20211125
                write(7,*) 'kswkw=',kswkw,'kresl=',kresl,'kmass=',kmass
                if(kswkw/=0) call propty
                if(kresl/=0) call stiff_u

                call dateandtime(curtime)   !20200220
                write(chkunit,2000)'Begin forming various matrix       at ',curtime  !20200220
                write(*,2000)      'Begin forming various matrix       at ',curtime   !20200220



                if(kmass/=0) call mcmatrx('U')
                if(ksmat/=0) call mcmatrx('W')
                if(allocated(fmass))call fmass_assemble

                if(khmat/=0) call hmatrx('W')
                if(kthmat/=0)call htmatrx

                print *,'ktsmat=',ktsmat
                if(ktsmat/=0)call stmatrx

                if(kthmat/=0)call assemble_boundt_estif
                if(kthmat/=0.and.algo_pipe==3)call assemble_pipe_estif
                if(kqmat/=0) call upwcouple
                if(stabpw==1.and.khmat/=0) call stabpatch   !!stablize

                if(kgrav/=0)call gravity
                if(kldfl/=0)call loadfl
                if (kresl/=0.or.ksmat/=0.or.khmat/=0.or.kqmat/=0.or.kmass/=0  &
                    .or.kthmat/=0.or.ktsmat/=0)then

                    call dateandtime(curtime)  !20200220
                    write(chkunit,2000)'Begin assemble golbal matrix       at ',curtime !20200220
                    write(*,2000)      'Begin assemble golbal matrix       at ',curtime !20200220


                    if(type_solver/='JPCG')global_stiff1=0.0
                    if(nonsym/=0.and.type_solver=='PROFILE')global_stiff2=0.0
                    if (type_solver=='JPCG')then
                        do ielem=1,nelem
                            element(ielem)%estif=0.0
                        end do
                    endif
                    call estif_assemble
                    call COUPLE_ASSEMBLE
                    if(nmcon/=0) call concentrated_mass_matrix !2013/4/11
                    if(stabpw==1) call stabpw_assemble           !! stablize
                    if(nifsgroup/=0 .and.(type_solver=='PROFILE'.or.type_solver=='PARDISO') )call assemble_interface_fluid_solid          !!ifs2000
                    if(nifsgroup/=0 .and.type_solver=='SSORPBCG')call assemble_interface_fluid_solid_SSORPBCG !!ifs2000
                    if(nabsfgroup/=0.and.(type_solver=='PROFILE'.or.type_solver=='PARDISO'))call assemble_absorb_fluid                   !!ifs2000
                    if(nabssgroup/=0.and.(type_solver=='PROFILE'.or.type_solver=='PARDISO'))call assemble_absorb_solid                   !!ifs2000
                    if(nabssgroup/=0.and.type_solver=='SSORPBCG')call assemble_absorb_solid_SSORPBCG          !!ifs2000
                    if((Icaddmass>=2.or.ifsnedge/=0).and.(type_solver=='PROFILE'.or.type_solver=='PARDISO') )call assemble_stiff_ifs2006   !!20220330 Li
                    if(ifsnedge/=0  .and.type_solver=='SSORPBCG')call assemble_stiff_ifs2006_SSORPBCG         !!ifs2006 zhao, 06/03/29

                    !do itotv=1,ntotv   !20221014
                    !   if (totveq(itotv)/=0)then
                    ! if(abs(global_stiff1(iseq(totveq(itotv)))).le.1.e-10.and.ifsnedge==0.and.nifsgroup==0.AND.TYPE_PROBLEM/='S')global_stiff1(iseq(totveq(itotv)))=1.e30
                    !   endif
                    !enddo  !20221014

                endif

                ! write(7,*)'global_stiff1***'
                !  do itotv=1,ntotv
                !if(totveq(itotv)/=0)   &
                ! write(7,*)itotv,totveq(itotv),global_stiff1(iseq(totveq(itotv)))
                !  end do


                if(iiter==1.or.kgrav/=0.or.kldfl/=0) call force_external
                !   write(chkunit,*)'itotv,tofor'
                !do itotv=1,ntotv
                !	if(abs(totveq(itotv)).ne.0) then
                !		write(chkunit,*)itotv,tofor(itotv)
                !	endif
                !end do


                if (iiter==1) then     ! iiter==1

                    call heat_flow_charge(0)          !!!20200316 pipe
                    call eload_initialize
                    call predict
                    call residu_f

                    !write(7,*)'eload1=',element(1)%field(1)%eload
                    !write(7,*)'eload2=',element(2)%field(1)%eload

                    call eload_couple
                    call eload_field  !20210417
                    if(algo_pipe>3) call pipe_cool_eload
                    call eload_interface_fluid_solid  !!ifs2000
                    call eload_absorb_fluid           !!ifs2000
                    call eload_absorb_solid           !!ifs2000
                    call eload_ifs2006                !!ifs2006 zhao, 06/03/29
                    if(stabpw==1)call stabload
                    call force_internal
                endif ! iiter=1    !should change for nssoil

                if(ngaps/=0.and.iblks>=iblks_bt.and.iiter==1)call ctfor_to_tofor(tofor0,tofor)  !!ctt2005

                if(type_nl==8.and.(iiter>1.or.(kstat==2.and.iiter>2))) goto 2  !MNR
                if ((type_solver=='PROFILE'.or.type_solver=='PARDISO').and.                     &
                    (kresl/=0.or.ksmat/=0.or.khmat/=0.or.kqmat/=0.or.kmass/=0  &
                    .or.kthmat/=0.or.ktsmat/=0)) then
                    call dateandtime(curtime)  !20200220
                    write(chkunit,2000)'Begin factorize global matrix      at ',curtime !20200220
                    write(*,2000)      'Begin factorize global matrix      at ',curtime	!20200220

                    operation='FACTORIZE'
                    call solve
                    !20230216 形成反演边界温度需要的C矩阵
                    if(nbackdT==2.and.istep==1) then
                        if(allocated(cmatrix_dtv))deallocate(cmatrix_dtv)
                        if(allocated(inv_cmatrix_dtv2))deallocate(inv_cmatrix_dtv2)
                        allocate(cmatrix_dtv(Npoints_pbx,nfixsets),inv_cmatrix_dtv2(nfixsets,nfixsets))
                        call cmatrix_dtv_formation(cmatrix_dtv,inv_cmatrix_dtv2)
                        !print *,'nfixsets=',nfixsets
                        !print *,'cmatrix_dtv=',cmatrix_dtv(:,1)
                        !print *,'inv_cmatrix_dtv2=',inv_cmatrix_dtv2

                    endif
                    !20230216
                end if
2               continue
                if (ngaps/=0.and.(iiter==1.and.istep==inc_step).and.iblks==iblks_bt)then   !!ctt2005
                    if (restart_ctt==0)then !restart_ctt
                        call forAdirect !fzx !形成A矩阵
                        kdimn=ndimn
                        if(block_stab==1)kdimn=3*(ndimn-1) !2015/11/17
                        allocate(rot(kdimn,kdimn))


                        do igapb=1,ngapb

                            if(block_appear_process(igapb,iblks)==0)cycle  !20200331
                            if(gapb(igapb)%eblock==0) cycle  !2017/11/19
                            npgblock=gapb(igapb)%npgblock

                            allocate(unitl(kdimn),unitg(kdimn))
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
                                    if(ndimn==3)rot(4:6,4:6)= rot(1:ndimn,1:ndimn)
                                endif
                                do idimn=1,kdimn
                                    itotvbt=(ipoin-1)*kdimn+idimn
                                    unitl=0.
                                    unitl(idimn)=1.*coef
                                    unitg=transpose(rot).x.unitl
                                    rvector=0.

                                    call  unit_force_trans(igapb,kdimn,ij,unitg,igaps,ipair,rvector)
                                    !
                                    operation='SOLVE'
                                    call solve

                                    do kpoin=1,npgblock
                                        jgaps=gapb(igapb)%nodegblock_igaps(kpoin)   !因为上面有idimn循环，这里不能用igaps变量
                                        jpair=gapb(igapb)%nodegblock_ipairs(kpoin)
                                        ij0=gapb(igapb)%nodegblock_onetwo(kpoin)    !因为上面有idimn循环，这里不能用ij变量      !！2017/04/03
                                        call result_node_to_center(kdimn,ij0,jgaps,jpair,result,unitg)

                                        do jdimn=1,kdimn
                                            jtotvbt=(kpoin-1)*kdimn+jdimn
                                            gapb(igapb)%cmatrix(jtotvbt,itotvbt)=unitg(jdimn)
                                        end do
                                    end do  !kpoin


                                    if(gapb(igapb)%nrdof>0)call dfat_rigid_ctfor(igapb,itotvbt,result)   !dfat due to ctfor

                                end do    ! idimn
                            end do    ! ipoin

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
                                    if(ndimn==3)rot(4:6,4:6)= rot(1:ndimn,1:ndimn)
                                endif
                                rot=coef*rot
                                allocate(cmatrixl(kdimn,kdimn))
                                do jpoin=1,npgblock
                                    jtotv=(jpoin-1)*kdimn
                                    itotv=(ipoin-1)*kdimn
                                    cmatrixl=gapb(igapb)%cmatrix(itotv+1:itotv+kdimn,jtotv+1:jtotv+kdimn)
                                    gapb(igapb)%cmatrix(itotv+1:itotv+kdimn,jtotv+1:jtotv+kdimn)=  &
                                        (rot.x.cmatrixl)
                                end do
                                deallocate(cmatrixl)
                            end do
                            !!!!!!!!!!!!!!!!!!!!!!!
                            !if(nonsbt/=0)then
                            do idimn=1,gapb(igapb)%nrdof
                                itotvbt=npgblock*kdimn+idimn
                                if(gapb(igapb)%listrdof(idimn)==0)cycle
                                allocate(force_rigid(ntotv))
                                call force_unit_rigid_accs(igapb,idimn,force_rigid)

                                operation='SOLVE'
                                call solve

                                do ipoin=1,npgblock
                                    !kpoin=gapb(igapb)%nodegblock(ipoin)
                                    igaps=gapb(igapb)%nodegblock_igaps(ipoin)
                                    ipair=gapb(igapb)%nodegblock_ipairs(ipoin)
                                    ij0=gapb(igapb)%nodegblock_onetwo(ipoin)   !2017/04/03
                                    coef=1.
                                    if(ij==2)coef=-1.
                                    rot=0.
                                    rot(1:ndimn,1:ndimn)=gaps(igaps)%rot(:,:,ipair)
                                    if(kdimn>ndimn)then
                                        if(ndimn==2)rot(3,3)=1.
                                        if(ndimn==3)rot(4:6,4:6)= rot(1:ndimn,1:ndimn)
                                    endif
                                    rot=coef*rot

                                    call result_node_to_center(kdimn,ij0,igaps,ipair,result,unitg)

                                    unitl=rot.x.unitg
                                    jtotv=(ipoin-1)*kdimn
                                    gapb(igapb)%cmatrix(jtotv+1:jtotv+kdimn,itotvbt)=   &
                                        gapb(igapb)%cmatrix(jtotv+1:jtotv+kdimn,itotvbt)+unitl
                                end do

                                call dfat_rigid(igapb,idimn,result,force_rigid)   !dfat due to stfor_inc
                                deallocate(force_rigid)
                            end do
                            !endif
                            deallocate(unitl,unitg)   !20121216
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!
                        end do  !igapb
                        do igapb=1,ngapb
                            if(block_appear_process(igapb,iblks)==0)cycle  !20200331
                            do itotvbt=1,gapb(igapb)%ntotv_bt
                                do jtotvbt=1,gapb(igapb)%ntotv_bt
                                    write(recttunit)gapb(igapb)%cmatrix(itotvbt,jtotvbt)
                                enddo
                            enddo
                        enddo  !igapb

                        deallocate(rot)
                    elseif(restart_ctt==1)then !restart_ctt
                        write(7,*)'read_cmatrix'
                        call forAdirect !fzx !形成A矩阵
                        rewind(recttunit)
                        do igapb=1,ngapb
                            if(block_appear_process(igapb,iblks)==0)cycle  !20200331
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

                if (nwcpipe/=0.and.(iiter==1.and.istep==inc_step))  &   !!20210411
                    call Tcmatrix_c_formation



90              format(10e12.5)

                !write(7,*)' itotv,tofor,stfor'
                !do itotv=1,ntotv
                !    write(7,*)itotv,tofor(itotv),stfor(itotv)
                !end do


                rvector=0.0
                !write(7,*)' itotv,totveq(itotv),rvector(totveq(itotv)),global_stiff1(iseq(totveq(itotv)))'

                if (type_solver/='JPCG') then
                    do itotv=1,ntotv
                        if (totveq(itotv)/=0)then
                            if (allocated(fexta))then  !!nstoks
                                rvector(totveq(itotv))=rvector(totveq(itotv))+tofor(itotv)+fexta(itotv)-stfor(itotv)
                            else
                                rvector(totveq(itotv))=rvector(totveq(itotv))+tofor(itotv)-stfor(itotv)

                                !if(abs(rvector(totveq(itotv)))>1.e-3) &
                                !write(7,*)itotv,totveq(itotv),rvector(totveq(itotv)),global_stiff1(iseq(totveq(itotv)))
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
                                if(iieq/=0)rvector(iieq)=rvector(iieq)+(tofor(itotv)-stfor(itotv))*trans(itotv)%rintf(iintf)
                            enddo
                        endif
                        !      if(totveq(itotv)/=0)then
                        !!if(abs(rvector(totveq(itotv)))>1.e-3) then
                        !!if(abs(totveq(itotv)).ne.0) then
                        !	!write(chkunit,'(2i8,20e16.6)')itotv,totveq(itotv),rvector(totveq(itotv)),global_stiff1(iseq(totveq(itotv)))
                        !      !endif
                        !      endif


                    enddo
                    !!int2000
                else !!if (type_solver/='JPCG') then

                    if (allocated(fexta))then   !!nstoks
                        rvector=tofor+fexta-stfor
                    else
                        rvector=tofor-stfor
                    endif
                endif


                if (type_nl==8) then
                    if(kstat/=2)call bfgsr(iiter)
                    if(kstat==2)call bfgsr(iiter-1)
                else
                    call dateandtime(curtime) !20200220
                    write(chkunit,2000)'Begin solution global equation     at ',curtime !20200220
                    write(*,2000)      'Begin solution global equation     at ',curtime !20200220
                    operation='SOLVE'
                    call solve
                endif
                !write(chkunit,'(a)')'result='
                !do itotv=1,ntotv
                !write(chkunit,*)itotv,result(itotv)
                !end do




                if(ngaps/=0.and.iblks>=iblks_bt)call solve_ctt  !!ctt2005
                if(nwcpipe/=0)call solve_heat_quantity_of_wc  !!20210411

                call varupdate !20230216

                !20230216 反演边界温度增量速率
                if(nbackdT==2.and.iiter==1) then
                    allocate(observstar(Npoints_pbx),dtv(nfixsets),dtvi(nfixsets))
                    observstar=0.

                    do i=1,mvalue
                        !print *,i,'Value_observ(i)%istep=',Value_observ(i)%istep
                        if(Value_observ(i)%iblks/=iblks)cycle
                        if(Value_observ(i)%iincs/=iincs)cycle
                        if(Value_observ(i)%istep/=istep)cycle
                        ivalue_point=Value_observ(i)%ivalue_point
                        observstar(ivalue_point)=Value_observ(i)%value_measure
                        nintf=para_points(ivalue_point)%nintf
                        listf=>para_points(ivalue_point)%listf
                        rintf=>para_points(ivalue_point)%rintf
                        Tpredict=dot_product(result_zero(listf),rintf)
                        !print *,'Tpredict=',Tpredict,'observstar(ivalue_point)=',observstar(ivalue_point)
                        observstar(ivalue_point)=observstar(ivalue_point)-Tpredict
                        nullify(listf,rintf)
                    end do
                    !print *,'observstar=',observstar
                    dTv=transpose(cmatrix_dtv).x.observstar
                    dtvi=inv_cmatrix_dtv2.x.dtv
                    !print *,'dtvi=',dtvi
                    dtvi=dtvi/(theta1*ditime)
                    !print *,'dtvi0=',dtvi


                    do i=1,ndofix    !20231130
                        mfixset=prescrib(i)%mfixset
                        do j=1,mfixset !20231130
                            ifixset=prescrib(i)%mlist(j)
                            idofn=prescrib(i)%ldofix
                            result_first(idofn)= result_first(idofn)+dtvi(ifixset)*prescrib(i)%rintf(j)
                            result_zero(idofn)= result_zero(idofn)+dtvi(ifixset)*prescrib(i)%rintf(j)*ditime !(theta1*ditime)
                        end do !20231130
                    end do !20231130


                    !   do i=1,ndofix
                    !ifixset=prescrib(i)%ifixset
                    !idofn=prescrib(i)%ldofix
                    !result_first(idofn)= result_first(idofn)+dtvi(ifixset)
                    !result_zero(idofn)= result_zero(idofn)+dtvi(ifixset)*ditime !(theta1*ditime)
                    !  end do

                    deallocate(observstar,dtv,dtvi)

                endif
                !20230216

                !call varupdate !20230216
                !call relative_dis_watertight !20231007 止水 !20240305
                call eload_initialize
                call residu_f
                call eload_couple
                call eload_field  !20210417
                if(algo_pipe>3) call pipe_cool_eload
                call eload_interface_fluid_solid  !!ifs2000
                call eload_absorb_fluid           !!ifs2000
                call eload_absorb_solid           !!ifs2000
                call eload_ifs2006                !!ifs2006 zhao, 06/03/29
                call heat_flow_charge(1)          !!! 20200316 pipe
                if(stabpw==1)call stabload


                call reaction_prescribed

                call conver_load
                if(nchek==0)call conver_nodal_value

                if(nchek==0)exit

            end do   !! loop for iiter
            ! call acc_rigid
            if(istatec==0) &
                call state_and_stiff_2021

            call local_stress    !20210207
            call contact_state(1) !20210207

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



            if(mdofn>7)then
                if (ifsnedge==0.and.type_problem=='F'.and.order_time_mdofn(8)==1)then
                    do ipoin=1,npoin
                        itotv=nodfn(lmdofn(8),ipoin)
                        IF(ITOTV/=0)result_second(itotv)=0.
                    enddo
                endif
            endif
            if(allocated(fexta))call find_fexta !!nstoks

            call cvoid(1)

            ! special for temperature (only) field
            if (outintw/=0)then
                resultm=0.
                do ipoin=1,npoin
                    idofn=nodfn(1,ipoin)
                    if(idofn/=0)resultm(ipoin)=result_zero(idofn)-result0(idofn)
                end do
                do ilink=1,ntlink
                    node1=tlink(1,ilink)
                    node2=tlink(2,ilink)
                    resultm(node1)=resultm(node2)
                end do
            endif
            ! end of special

            call gpvarupdate
            !write(7,*)'af gpvarupdate','stres0=',element(1)%field(1)%gpvar(1:3,1) !,'stres=',element(1)%field(1)%gpvar(1:3,1)

            if(gamamax/=0)call gamamaxupdate !20231125YL 更新地震过程中最大动剪应变
            if (istep/noutn*noutn==istep)then
                iwriten=iwriten+1
                call out_record
                call outputres !for output
                call out_full_write
            endif

            if (istep/noutf*noutf==istep)then

                if(outplot(1:3)=='GID')   call OUT_GID_WRITE
                if(outplot(1:6)=='COSMOS')call OUT_COSMOS_WRITE
            endif

            if (outintw/=0) then
                !write(7,*)'resultm(931)=',resultm(931),'result_zero(931)=',result_zero(931)
                write(outint)resultm
            endif

            if(istep/nresta*nresta==istep)call resta_read_write(-1)

            if(nstepjq>0)then !20231215YL
                if(istep/nstepjq*nstepjq==istep)call liquifaction_judge
                if(istep/nstepJP*nstepJP==istep)call permdeform_judge
            endif !20231215YL

            if(nforce/=0.or.ngaps/=0)call force_interface
            !         if(nforce/=0.or.ngaps/=0)call write_force_interface
            if((nforce/=0.or.ngaps/=0).and.nextrf==0)call write_force_interface
            if(kstab/=0.)call safety_factor  !2019/03/20

            if (nextrf/=0)then  !2004/9/11
                do idofix=1,ndofix
                    itcurve =prescrib(idofix)%itcurve
                    type_curve=tcurves(itcurve)%type_curve
                    if (type_curve=='EXTRAPOLATION')then
                        if (istep>nextrf)then
                            do iextrf=1,nextrf-1
                                prescrib(idofix)%value_ext(:,iextrf+1)=prescrib(idofix)%value_ext(:,iextrf)
                            end do
                            prescrib(idofix)%value_ext(:,1)=result_zero(prescrib(idofix)%listep)
                        else
                            if(istep==1)prescrib(idofix)%value_ext(:,1)=0.
                            do iextrf=1,istep-1
                                prescrib(idofix)%value_ext(:,iextrf+1)=  &
                                    prescrib(idofix)%value_ext(:,iextrf)
                            end do
                            prescrib(idofix)%value_ext(:,1)=result_zero(prescrib(idofix)%listep)
                        endif
                    endif
                end do
            endif !2004/9/11
            !hxl2006 MIF
            if(type_ABC=='MIF')then
                do idofix=1,ndofix   !sanshe  hxl
                    ldofixb=>prescrib(idofix)%ldofixb
                    do ilaymif=1,nlaymif
                        !dissanru(istep,ldofixb(i))=result_zero(ldofixb(i))-inpru(ldofixb(i))
                        !dissanzi(istep,ldofixb(i))=result_zero(ldofixb(i))-inpzi(ldofixb(i))
                        disA(ldofixb(ilaymif))=result_zero(ldofixb(ilaymif))-inpru(ldofixb(ilaymif))
                        !对底边界，将总波场分解为入射波场与散射波场，disA即为散射波
                        disB(ldofixb(ilaymif))=result_zero(ldofixb(ilaymif))-inpzi(ldofixb(ilaymif))
                        !对侧边界，将总波场分解为自由波场与散射波场，disB即为散射波

                        !原因在于透射边界仅对散射波场，保证散射波场能够穿过人工边界而透向无限远处，但同时与要允许入射波场能够向上传播
                        !故进行波场分离。
                    end do
                    nullify(ldofixb)
                end do
            endif
            deallocate(inpru,inpzi)
            ! end hxl2006 MIF

            if(cdtest>=1)then
                nstre =group(1)%nstre
                call gpq (element(1)%field(1)%gpvar(1:nstre,1),p,q,eta)
                CVV=Q
                xl0=0.
                DIST0=ABS(CVV-XL0)
                DIST1=ABS(CVV-XL1(ITRAM))
                DIST2=ABS(CVV-XL2(ITRAM))
                DIST =ABS(XL1(ITRAM)-XL2(ITRAM))
                DISTF=ABS(XL2(ITRAM)-XL0)
                ICREV=0
                ICEND=0
                IF((DIST1.GE.DIST).OR.(DIST2.GE.DIST)) ICREV=1
                !C------------------------------------------------------
                !C       CHECK END OF TRAM
                !C-----------------------------------------------------
                NC1=2*NCYC(ITRAM)-1
                NC2=2*NCYC(ITRAM)
                !	    print *,'iseg=',iseg,'icrev=',icrev,'fincre=',fincre
                IF((cdtest==1).AND.(ISEG.EQ.NC1).AND.(ICREV.EQ.1)) ICEND=1
                IF((cdtest==2).AND.(ISEG.EQ.NC2).AND.(DIST2.GE.DISTF)) ICEND=1
                IF(ICEND.EQ.1) ITRAM=ITRAM+1
                IF((ICREV.NE.1).AND.(ICEND.NE.1)) GO TO 32
                ISEG=ISEG+1
                IF(ICEND.EQ.1)  ISEG=0
                fincre=-fincre
32              CONTINUE
                !       print *,'icend=',icend,'itram=',itram
                IF((ICEND.EQ.1).AND.(ITRAM.GT.NTRAM)) stop
            endif

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
                        !write(7,*)'ivalue=',ivalue,'idofn=',idofn,'nodfn(idofn,listf)=',nodfn(idofn,listf)
                        !write(7,*)'ivalue_point=',ivalue_point,'listf=',listf,'Value_observ(ivalue)%value_computation=',Value_observ(ivalue)%value_computation


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

            if(outinp<0)then  !非稳定渗流场分析时向oip文件输出结点压力 20220623
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

            if(upliftin<0)then  !非稳定渗流场分析时向upf文件输出结点压力 20220623
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

            endif
            !if((bparameter>=1.and.bparameter<=2).and.balgor>=1) call dudx !20230430

        end do       !! for istep

        deallocate(disA,disB,disA_1,disB_1,disA_2,disB_2) !hxl2006 MIF
    end do     !! loop for iincs

    if(cdtest>=1)deallocate(xl1,xl2,ncyc,ic)
    if(outintw/=0)deallocate(resultm,result0)
    if(ngaps/=0) deallocate(tofor0)
2000 format(a,a)
    end subroutine time_dependent

    subroutine value_submodel_boundary  !20210324

    integer(ink) iu,it,itotv,idimn,kkdimn
    real   (irk), allocatable::value(:,:)

    if(res_u==1)then

        kkdimn=ndimn
        if(res_rot/=0) kkdimn=3*(ndimn-1)  !20221202

        allocate(value(kkdimn,tbpointsu))
        read(resbunit)value


        do iu=1,tbpointsu
            ipoin=listbpointsu_t(iu)
            do idimn=1,kkdimn
                itotv=nodfn(idimn,ipoin)
                if(istep==1) then
                    result_zero(itotv)=value(idimn,iu)
                    deltafi(itotv)=value(idimn,iu)
                else
                    deltafi(itotv)=value(idimn,iu)-result_zero(itotv)
                    result_zero(itotv)=value(idimn,iu)
                endif
                !write(chkunit,*)'iu=',iu,'ipoin=',ipoin,'idimn=',idimn,'value=',deltafi(itotv)

            end do
        end do


        if(res_v==1)then
            read(resbunit)value
            do iu=1,tbpointsu
                ipoin=listbpointsu_t(iu)
                do idimn=1,kkdimn
                    itotv=nodfn(idimn,ipoin)
                    result_first(itotv)=value(idimn,iu)
                end do
            end do
        endif

        if(res_a==1)then
            read(resbunit)value
            do iu=1,tbpointsu
                ipoin=listbpointsu_t(iu)
                do idimn=1,kkdimn
                    itotv=nodfn(idimn,ipoin)
                    result_second(itotv)=value(idimn,iu)
                end do
            end do
        endif
        deallocate(value)
    endif


    if(res_T==1)then
        allocate(value(1,tbpointst))
        read(resbunit)value

        do it=1,tbpointsT
            ipoin=listbpointst_t(it)
            itotv=nodfn(lmdofn(10),ipoin)
            result_zero(itotv)=value(1,it)
        end do
        if(res_Tv==1)then
            read(resbunit)value

            do it=1,tbpointsT
                ipoin=listbpointst_t(it)
                itotv=nodfn(lmdofn(10),ipoin)
                result_first(itotv)=value(1,it)
            end do
        endif
        deallocate(value)
    endif

    !write(7,*)'pressure='
    if(res_P==1)then
        allocate(value(1,tbpointsp))
        read(resbunit)value
        do it=1,tbpointsP
            ipoin=listbpointsp_t(it)
            itotv=nodfn(lmdofn(8),ipoin)
            result_zero(itotv)=value(1,it)
            !write(7,*)ipoin,result_zero(itotv)
        end do
        !write(7,*)'pressure_v='

        if(res_Pv==1)then
            read(resbunit)value
            do it=1,tbpointsP
                ipoin=listbpointsp_t(it)
                itotv=nodfn(lmdofn(8),ipoin)
                result_first(itotv)=value(1,it)
                !write(7,*)ipoin,result_first(itotv)

            end do
        endif
        !write(7,*)'pressure_a='

        if(res_Pa==1)then
            read(resbunit)value
            do it=1,tbpointsP
                ipoin=listbpointsp_t(it)
                itotv=nodfn(lmdofn(8),ipoin)
                result_second(itotv)=value(1,it)
                !write(7,*)ipoin,result_second(itotv)

            end do
        endif
        deallocate(value)
    endif

    end subroutine value_submodel_boundary !20210324

    SUBROUTINE force_unit_rigid_accs(igapb,jdimn,force_rigid)
    character(10)fieldid
    integer(ink) igroup,index,nnode_f, nevab_f, ic,ielgroup, ielem,ievab,itotv,igapb,jdimn,ipoin,jpoin,jgapb,kdimn,kkdimn
    real   (irk)  coef,alfa,beta,force_rigid(:)
    real   (irk), allocatable::fstif(:,:),eload(:),value(:),result_second_rigid(:)
    real   (irk), pointer::fstif0(:,:)
    integer(ink), pointer::ldofs(:)

    kkdimn=ndimn
    if(block_stab==1)kkdimn=3*(ndimn-1) !2015/11/17
    allocate(result_second_rigid(ntotv))

    rvector=0.;force_rigid=0.;result_second_rigid=0.
    do jpoin=1,gapb(igapb)%npblock
        ipoin=gapb(igapb)%nodeblock(jpoin)
        do kdimn=1,kkdimn  !idimn
            itotv=nodfn(kdimn,ipoin)
            if(itotv/=0) &
                result_second_rigid(itotv)=1.*gapb(igapb)%npdisp(kdimn,jpoin,jdimn)
        end do   !kdimn
    end do   !jpoin

    !	write(7,*)'result_second_rigid=',result_second_rigid

    do jgapb=1,gapb(igapb)%ngroupb
        igroup=gapb(igapb)%listgroupb(jgapb)
        fieldid=group(igroup)%fieldid
        if(appear(igroup)>0)  then
            ! get information from the group level
            index    = group(igroup)%index
            alfa=group(igroup)%alfa
            beta=group(igroup)%beta
            nnode_f = elkn(index)%el_field(1)%nnode_f
            nevab_f = nnode_f*group(igroup)%dof(1)%nfdof
            allocate(fstif(nevab_f,nevab_f),eload(nevab_f),value(nevab_f))
            ! loop for k(h) and m(c)
            DO ielgroup = 1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                if(fieldid(1:1)=='U') then
                    if(associated(element(ielem)%field(1)%khandmc(2)%fstif)) then
                        coef=1.0+alfa*beeta1*ditime                         !-------------------------------!
                        ldofs=>element(ielem)%field(1)%ldofs_f
                        value=result_second_rigid(ldofs)
                        fstif0=>element(ielem)%field(1)%khandmc(2)%fstif
                        ic=size(fstif0,dim=2)
                        if(ic==1)then
                            fstif=0.0
                            do ievab=1,nevab_f
                                fstif(ievab,ievab)=fstif0(ievab,1)
                            end do
                        else
                            fstif=fstif0
                        end if
                        fstif=coef*fstif
                        eload=fstif.x.value
                        force_rigid(ldofs)=force_rigid(ldofs)+eload
                        nullify(fstif0,ldofs)
                    endif
                endif
            end do       !!ielgroup
            deallocate(fstif,eload,value)
        end if    !! for do while
    end do     !!  for igroup

    do itotv=1,ntotv
        if (totveq(itotv)>0)rvector(totveq(itotv))=rvector(totveq(itotv))-force_rigid(itotv)
    end do


    deallocate(result_second_rigid)

    END SUBROUTINE force_unit_rigid_accs

    SUBROUTINE dfat_rigid(igapb,idimn,resi,force_rigid)
    character(10)fieldid
    character(30)material
    integer(ink) igroup, nrfields, ifield,  index, order_time,     &
        nnode_f, nevab_f,   ic,   idimn,igapb, jgapb, npblock,npgblock,kdimn,       &
        ielgroup, ielem,   ievab,  ikh, anevab,  matno, nstre,jdimn,jtotvbt,itotvbt

    real   (irk)  coef,alfa,beta,lamda,resi(:),force_rigid(:)
    real   (irk), allocatable::fstif(:,:),eload(:),value(:),stfor_inc(:),df(:),dfat(:),mtrxA(:,:)
    real   (irk), pointer::fstif0(:,:)
    integer(ink), pointer::ldofs(:)

    allocate(stfor_inc(ntotv))
    stfor_inc=0.
    do jgapb=1,gapb(igapb)%ngroupb
        igroup=gapb(igapb)%listgroupb(jgapb)
        if(appear(igroup)>0)  then

            ! get information from the group level
            nrfields=group(igroup)%nrfields
            fieldid=group(igroup)%fieldid
            matno = group(igroup)%matno
            nstre=  group(igroup)%nstre
            if(fieldid(1:1)=='U') &
                material=props(matno)%mechanical%solid%material
            index    = group(igroup)%index
            alfa=group(igroup)%alfa
            beta=group(igroup)%beta
            do ifield=1,nrfields
                nnode_f = elkn(index)%el_field(ifield)%nnode_f
                nevab_f = nnode_f*group(igroup)%dof(ifield)%nfdof
                allocate(fstif(nevab_f,nevab_f),eload(nevab_f),value(nevab_f))
                ! loop for k(h) and m(c)
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    if(fieldid(ifield:ifield)=='U') then
                        do ikh=1,2
                            if(associated(element(ielem)%field(ifield)%khandmc(ikh)%fstif)) then
                                if(ikh==1)coef=beeta2*ditime**2+beta*beeta1*ditime          !
                                if(ikh==2)coef=1.0+alfa*beeta1*ditime                         !-------------------------------!

                                ldofs=>element(ielem)%field(ifield)%ldofs_f
                                value=resi(ldofs)
                                fstif0=>element(ielem)%field(ifield)%khandmc(ikh)%fstif
                                ic=size(fstif0,dim=2)
                                if(ikh==2.and.ic==1)then
                                    fstif=0.0
                                    do ievab=1,nevab_f
                                        fstif(ievab,ievab)=fstif0(ievab,1)
                                    end do
                                else
                                    fstif=fstif0
                                end if
                                fstif=coef*fstif
                                eload=fstif.x.value
                                stfor_inc(ldofs)=stfor_inc(ldofs)+eload
                                nullify(fstif0,ldofs)
                            endif
                        end do        !!end do ikh
                    endif
                end do       !!ielgroup
                deallocate(fstif,eload,value)
            end do     !! end do ifield
        end if    !! for do while
    end do     !!  for igroup

    kdimn=ndimn
    if(block_stab==1)kdimn=3*(ndimn-1) !2015/11/17
    allocate(df(kdimn),dfat((ndimn-1)*3),mtrxA(kdimn,(ndimn-1)*3))
    df=0. ;  dfat=0. ;mtrxA=0.

    dfat=0.
    npblock=gapb(igapb)%npblock
    npgblock=gapb(igapb)%npgblock
    do jpoin=1,gapb(igapb)%npblock
        ipoin=gapb(igapb)%nodeblock(jpoin)
        df=0.
        do jdimn=1,kdimn !严格的说，应该用cdofn,ndof
            itotv=nodfn(jdimn,ipoin)
            if(itotv/=0) &
                df(jdimn)=-stfor_inc(itotv)-force_rigid(itotv)
        enddo

        mtrxa=0.
        do jdimn=1,gapb(igapb)%nrdof
            mtrxa(:,jdimn)=gapb(igapb)%npdisp(:,jpoin,jdimn)
        end do

        dfat=dfat+matmul(transpose(mtrxA),df)

    enddo !jpoin
    jtotvbt=npgblock*kdimn+idimn
    do jdimn=1,gapb(igapb)%nrdof
        itotvbt=npgblock*kdimn+jdimn
        gapb(igapb)%cmatrix(itotvbt,jtotvbt)=gapb(igapb)%cmatrix(itotvbt,jtotvbt)+dfat(jdimn)

    enddo !jdimn

    deallocate(stfor_inc,df,dfat,mtrxA)

    END SUBROUTINE dfat_rigid

    SUBROUTINE dfat_rigid_ctfor(igapb,itotvbt,resi)
    character(10)fieldid
    character(30)material
    integer(ink) kpoin,ipoin,igroup, nrfields, ifield,  index, order_time,     &
        nnode_f, nevab_f,   ic,  igapb, jgapb, npblock,npgblock,kdimn,       &
        ielgroup, ielem,   ievab,  ikh, anevab,  matno, nstre,jdimn,jtotvbt,itotvbt

    real   (irk)  coef,alfa,beta,lamda,resi(:)
    real   (irk), allocatable::fstif(:,:),eload(:),value(:),stfor_inc(:),df(:),dfat(:),mtrxA(:,:)
    real   (irk), pointer::fstif0(:,:)
    integer(ink), pointer::ldofs(:)

    allocate(stfor_inc(ntotv))
    stfor_inc=0.

    do jgapb=1,gapb(igapb)%ngroupb
        igroup=gapb(igapb)%listgroupb(jgapb)
        if(appear(igroup)>0)  then

            ! get information from the group level
            nrfields=group(igroup)%nrfields
            fieldid=group(igroup)%fieldid
            matno = group(igroup)%matno
            nstre=  group(igroup)%nstre
            if(fieldid(1:1)=='U') &
                material=props(matno)%mechanical%solid%material
            index    = group(igroup)%index
            alfa=group(igroup)%alfa
            beta=group(igroup)%beta
            do ifield=1,nrfields
                nnode_f = elkn(index)%el_field(ifield)%nnode_f
                nevab_f = nnode_f*group(igroup)%dof(ifield)%nfdof
                allocate(fstif(nevab_f,nevab_f),eload(nevab_f),value(nevab_f))
                ! loop for k(h) and m(c)
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    if(fieldid(ifield:ifield)=='U') then
                        do ikh=1,2
                            if(associated(element(ielem)%field(ifield)%khandmc(ikh)%fstif)) then
                                if(ikh==1)coef=beeta2*ditime**2+beta*beeta1*ditime          !
                                if(ikh==2)coef=1.0+alfa*beeta1*ditime                         !-------------------------------!

                                ldofs=>element(ielem)%field(ifield)%ldofs_f
                                value=resi(ldofs)
                                fstif0=>element(ielem)%field(ifield)%khandmc(ikh)%fstif
                                ic=size(fstif0,dim=2)
                                if(ikh==2.and.ic==1)then
                                    fstif=0.0
                                    do ievab=1,nevab_f
                                        fstif(ievab,ievab)=fstif0(ievab,1)
                                    end do
                                else
                                    fstif=fstif0
                                end if
                                fstif=coef*fstif
                                eload=fstif.x.value
                                stfor_inc(ldofs)=stfor_inc(ldofs)-eload
                                nullify(fstif0,ldofs)
                            endif
                        end do        !!end do ikh
                    endif
                end do       !!ielgroup
                deallocate(fstif,eload,value)
            end do     !! end do ifield
        end if    !! for do while
    end do     !!  for igroup

    kdimn=ndimn
    if(block_stab==1)kdimn=(ndimn-1)*3

    allocate(df(kdimn),dfat((ndimn-1)*3),mtrxA(kdimn,(ndimn-1)*3))
    df=0. ; dfat=0. ; mtrxA=0.

    dfat=0.
    npblock=gapb(igapb)%npblock
    npgblock=gapb(igapb)%npgblock
    do jpoin=1,gapb(igapb)%npblock
        ipoin=gapb(igapb)%nodeblock(jpoin)
        df=0.
        do jdimn=1,kdimn !严格的说，应该用cdofn,ndof
            itotv=nodfn(jdimn,ipoin)
            if(itotv/=0) &
                df(jdimn)=stfor_inc(itotv)
        enddo

        mtrxa=0.
        do jdimn=1,gapb(igapb)%nrdof
            mtrxa(:,jdimn)=gapb(igapb)%npdisp(:,jpoin,jdimn)
        end do

        dfat=dfat+matmul(transpose(mtrxA),df)
    enddo !jpoin

    do jdimn=1,gapb(igapb)%nrdof
        jtotvbt=npgblock*kdimn+jdimn
        gapb(igapb)%cmatrix(jtotvbt,itotvbt)=gapb(igapb)%cmatrix(jtotvbt,itotvbt)+dfat(jdimn)
    enddo !jdimn

    deallocate(stfor_inc,df,dfat,mtrxA)

    END SUBROUTINE dfat_rigid_ctfor

    subroutine explicit

    character(80)text
    integer(ink) itotv,ielem,iintf,iieq
    integer(ink),allocatable::earthquake_curve(:)
    real   (irk), allocatable::rmid(:)
    real   (irk) time
    real   (irk),pointer::ymass(:,:)
    integer(ink),pointer::ldofs(:)
    if(allocated(rmid))deallocate(rmid)
    if(allocated(result))deallocate(result)
    if(allocated(rvector))deallocate(rvector)
    allocate(rmid(ntotv),result(ntotv),rvector(ntotv))

    print *,     'in explicit'
    if(iblks==1)call cvoid(0)
    allocate(earthquake_curve(ndimn),fachv(ndimn))
    read(mainunit,*)text
    read(mainunit,*)nincs,earthquake_curve(1:ndimn)
    time=0.0
    do iincs=1,lincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta
        read(mainunit,*)toler_force,toler_var(1:mdofn)
    end do

    do iincs=lincs+1,nincs
        print *,'  time depend     iincs=',  iincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta
        read(mainunit,*)toler_force,toler_var(1:mdofn)
        do istep=1,nstep
            print *, '     istep=',istep
            if(outintr>0.and.iblks>=outintr)trstep=trstep+1  !20200226
            time=time+ditime
            ttime=ttime+ditime        !! only for output

            call dfact_time_curve(ttime)

            where(earthquake_curve==0)
                fachv=0.0
            elsewhere
                fachv=tcurves(earthquake_curve)%dfact
            endwhere

            call modf_var_prescribed
            call heat_internal        ! temperature
            deltafi=0.0

            do iiter=1,miter

                call algort
                call porepr
                if(kswkw/=0)call propty
                if(kresl/=0)call stiff_u
                if(kmass/=0)call mcmatrx('U')
                if(ksmat/=0)call mcmatrx('W')
                if(ktsmat/=0)call stmatrx
                if(khmat/=0)call hmatrx('W')
                if(kthmat/=0)call htmatrx
                call assemble_boundt_estif
                if(kqmat/=0)call upwcouple
                if(stabpw==1.and.khmat/=0) call stabpatch   !!stablize
                !!explicit
                print *,'kmass=',kmass
                if (ktsmat/=0.or.kmass/=0) then
                    rmid=0.
                    do ielem=1,nelem
                        ymass=>element(ielem)%field(1)%khandmc(2)%fstif
                        ldofs=>element(ielem)%field(1)%ldofs_f
                        rmid(ldofs)=rmid(ldofs)+ymass(:,1)
                    end do
                end if

                !!int2000
                do itotv=1,ntotv
                    nintf=trans(itotv)%nintf
                    if (nintf/=0) then
                        do iintf=1,nintf
                            iieq=trans(itotv)%listf(iintf)
                            if(iieq/=0) &
                                rmid(iieq)=rmid(iieq)+rmid(itotv)*trans(itotv)%rintf(iintf)
                        end do
                    endif
                end do
                !!int2000

                !! end explicit

                if(kgrav/=0)call gravity
                if(kldfl/=0)call loadfl

                if(iiter==1) call force_external

                if (iiter==1) then     ! iiter==1

                    call eload_initialize
                    call predict
                    call residu_f
                    call eload_couple
                    call eload_field
                    if(stabpw==1)call stabload
                    call force_internal

                endif          ! iiter=1

                rvector=tofor-stfor

                !!int2000
                do itotv=1,ntotv
                    nintf=trans(itotv)%nintf
                    if (nintf/=0) then
                        do iintf=1,nintf
                            iieq=trans(itotv)%listf(iintf)
                            if(iieq/=0) &
                                rvector(iieq)=rvector(iieq)+(tofor(itotv)-stfor(itotv))*trans(itotv)%rintf(iintf)
                        end do
                    endif
                end do
                do itotv=1,ntotv
                    nintf=trans(itotv)%nintf
                    if(nintf==0) &
                        result(itotv)=rvector(itotv)/rmid(itotv)  !!!
                end do

                !!int2000
                do itotv=1,ntotv
                    nintf=trans(itotv)%nintf
                    if (nintf/=0) then
                        result(itotv)=0.
                        do iintf=1,nintf
                            iieq=trans(itotv)%listf(iintf)
                            result(itotv)=result(itotv)+result(iieq)*trans(itotv)%rintf(iintf)
                        enddo
                    endif
                end do
                !!int2000
                !!int2000
                !            do itotv=1,ntotv
                !               result(itotv)=rvector(itotv)/rmid(itotv)  !!!
                !            end do

                call varupdate

                call eload_initialize
                call residu_f
                call eload_couple
                call eload_field
                if(stabpw==1)call stabload
                call reaction_prescribed
                call conver_load
                if(nchek==0)call conver_nodal_value
                if(nchek==0)exit

            end do   !! loop for iiter

            call cvoid(1)
            call gpvarupdate

            if (istep/noutn*noutn==istep)then
                iwriten=iwriten+1
                call out_record
            endif

            if (istep/noutf*noutf==istep)then
                call out_full_write
                if(outplot(1:3)=='GID')   call OUT_GID_WRITE
                if(outplot(1:6)=='COSMOS')call OUT_COSMOS_WRITE
            endif
            if(istep/nresta*nresta==istep)call resta_read_write(-1)
        end do       !! for istep
    end do     !! loop for iincs

    end subroutine explicit
