    subroutine frequency_analysis !freq2006

    character(80)text,type_curve
    integer(ink) itotv,ielem
    integer(ink) ilink,node1,node2,ipoin,idofn
    real   (irk), allocatable::resultm(:)
    real   (irk) time
    integer(ink) iintf,nintf,iieq   !!int2000


    if(allocated(fachv))deallocate(fachv)
    allocate(fachv(ndimn))
    read(mainunit,*)text
    read(mainunit,*)nincs,fachv

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
            print *, 'iblks=',iblks,'iincs=',iincs,'istep=',istep
            time=time+ditime
            ttime=ttime+ditime        !! only for output
            print *,'ttime=',ttime

            !deltafi=0.0
            stforw=0.
            call modf_var_prescribed_w
            do iiter=1,miter
                print *,'iiter=',iiter
                if (istep==1.and.iiter==1)then
                    call stiff_u
                    call mcmatrx('U')
                    call mcmatrx('W')
                    if(allocated(fmass))call fmass_assemble
                    call hmatrx('W')
                    call upwcouple
                    !               if(stabpw==1) call stabpatch    !!stablize
                endif

                if (iiter==1)then
                    global_stiff1w=0.0
                    if(nonsym/=0)global_stiff2w=0.0
                    call estif_assemble_w
                    call couple_assemble_w
                    if(stabpw==1)     call stabpw_assemble_w           !! stablize
                    if(nifsgroup/=0)  call assemble_interface_fs_w !!ifs2000
                    if(nabsfgroup/=0) call assemble_absorb_fluid_w !!ifs2000
                    if(nabssgroup/=0) call assemble_absorb_solid_w !!ifs2000
                    if(ifsnedge/=0)   call assemble_stiff_ifs2006_w !ifs2006
                endif

                if (iiter==1) call force_external_w
                if (iiter==1)then
                    operation='FACTORIZE'
                    call solve
                end if

                rvectorw=(0.0,0.)
                do itotv=1,ntotv
                    if(totveq(itotv)/=0)rvectorw(totveq(itotv))=rvectorw(totveq(itotv))+ &
                        (toforw(itotv)-stforw(itotv))
                end do

                !!int2000
                do itotv=1,ntotv
                    nintf=trans(itotv)%nintf
                    if (nintf/=0) then
                        iieq=totveq(itotv)
                        if(iieq/=0)rvectorw(iieq)=(0.,0.)
                        do iintf=1,nintf
                            iieq=totveq(trans(itotv)%listf(iintf))
                            if(iieq/=0)rvectorw(iieq)=rvectorw(iieq)+(toforw(itotv)-stforw(itotv))*trans(itotv)%rintf(iintf)
                        end do
                    endif
                end do
                !!int2000
                !            write(7,*)'global_stif,rvector'
                !            do itotv=1,ntotv
                !               if (abs(totveq(itotv)).ne.0) then
                !                  write(chkunit,*)itotv,totveq(itotv),global_stiff1w(iseq(totveq(itotv))),rvectorw(totveq(itotv))
                !               endif
                !            end do



                operation='SOLVE'
                call solve

                stforw=(0.,0.)
                call varupdate_w
                call eload_couple_w
                call eload_field_w
                call eload_interface_fs_w      !!ifs2000
                call eload_absorb_fluid_w      !!ifs2000
                call eload_absorb_solid_w      !!ifs2000
                call eload_ifs2006_w           !ifs2006
                call conver_load_w
                if(nchek==0)exit

            end do   !! loop for iiter


            !         call gpvarupdate

            if (istep/noutn*noutn==istep)then
                iwriten=iwriten+1
                call out_record
            endif

            !         if(istep/noutf*noutf==istep)then
            !         call out_full_write
            if(outplot(1:3)=='GID')   call out_gid_write_w
            !         if(outplot(1:6)=='COSMOS')call out_cosmos_write
            !         endif

            if(istep/nresta*nresta==istep)call resta_read_write(-1)

            if(nforce/=0.or.ngaps/=0)call force_interface
            !         if(nforce/=0.or.ngaps/=0)call write_force_interface

        end do       !! for istep
    end do     !! loop for iincs

    end subroutine frequency_analysis

    subroutine response_spectrum ! only available for u field + w field

    character(80)text
    integer(ink) itotv,iieq,nintf,iintf,jtotv,idimn,imcon,  &
        jstep,max_dofn,max_poin,max_dim,jiter,ipoin,idofn,maloc(1)
    integer(ink),pointer::listf(:)
    real   (irk), allocatable::rmid(:),alfa(:),beta(:),lamda(:),deltafp(:),   &
        dis_shape(:,:),mdelta(:),disone(:),loadone(:), &
        disthree(:),loadthree(:),kdelta(:),   &
        distwo(:),loadtwo(:) !zhao
    real   (irk),pointer::rintf(:)
    real   (irk) time,lamda_iter,omega,alfa_mid,beta_mid,coefx,coefy,coefz,kstar,coef !zhao


    print *,'response_spectrum'
    if(allocated(rmid))deallocate(rmid)
    if(allocated(rvector))deallocate(rvector)
    allocate(rmid(ntotv),rvector(ntotv),mdelta(ntotv),disone(ntotv),loadone(ntotv),kdelta(ntotv))
    if(nifsgroup/=0)allocate(deltafp(ntotv))
    allocate(disthree(ntotv),loadthree(ntotv),distwo(ntotv),loadtwo(ntotv))

    read(mainunit,*)text
    read(mainunit,*)nincs
    time=0.0
    do iincs=1,lincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta,nmcon !20210118
        read(mainunit,*)toler_force,toler_var(1:mdofn)
    end do

    do iincs=lincs+1,nincs
        print *,'  time depend     iincs=',  iincs
        read(mainunit,*)miter,ditime,noutn,noutf,nstep,inc_step,nresta,nmcon  !20210118
        read(mainunit,*)toler_force,toler_var(1:mdofn)
        read(mainunit,*)max_poin,max_dim
        if(nmcon/=0)allocate(lmcon(nmcon),rmcon(ndimn,nmcon))
        if (nmcon/=0) then
            read(mainunit,*)coef
            do imcon=1,nmcon
                read(mainunit,*)i0,lmcon(imcon),rmcon(:,imcon)
            end do
            rmcon=abs(rmcon)*coef
        endif
        max_dofn=nodfn(max_dim,max_poin)

        allocate(beta(nstep),alfa(nstep),dis_shape(ntotv,nstep),lamda(nstep))

        do istep=1,nstep

            print *, 'istep=',istep
            if(outintr>0.and.iblks>=outintr)trstep=trstep+1  !20200226
            time=time+ditime
            ttime=ttime+ditime        !! only for output
            call dfact_time_curve(ttime)
            call modf_var_prescribed

            if (istep==1) then
                iiter=1
                call stiff_u
                call mcmatrx('U')
                call hmatrx('W')
                global_stiff1=0.
                call estif_assem_response
                !write(7,*)'global_stiff1**'
                !do itotv=1,neq
                !    write(7,*)itotv,global_stiff1(iseq(itotv))
                !   if(abs(global_stiff1(iseq(itotv))).le.1.e-5)global_stiff1(iseq(itotv))=1.e30
                !end do
                !            call skfacs_layer(global_stiff1,iseq,0)
                if (type_solver=='PROFILE') then !zhao 20070829
                    call skfacs_layer(global_stiff1,iseq,1)
                    call skfacs_layer(global_stiff1,iseq,2)
                    write(7,*)'neq_layer1=',neq_layer1,'neq=',neq
                ELSEIF(type_solver=='PARDISO')THEN
                    operation='FACTORIZE'
                    call solve
                endif
            endif

            if (istep==1) then
                deltafi=0.
                disone=0.
                distwo=0. !zhao
                disthree=0.

                do ipoin=1,npoin
                    DO idimn=1,ndimn
                        itotv=nodfn(idimn,ipoin)
                        if (itotv/=0) then
                            if(iffix(itotv)==0)deltafi(itotv)=1.
                        endif
                    enddo
                end do

                do ipoin=1,npoin
                    itotv=nodfn(1,ipoin)
                    if (itotv/=0) then
                        if(iffix(itotv)==0)disone(itotv)=1.
                    endif
                end do

                if (cdofn>1)then
                    do ipoin=1,npoin
                        itotv=nodfn(2,ipoin)
                        if (itotv/=0) then
                            if(iffix(itotv)==0)distwo(itotv)=1.
                        endif
                    end do
                endif

                if (cdofn>2)then
                    do ipoin=1,npoin
                        itotv=nodfn(3,ipoin)
                        if (itotv/=0) then
                            if(iffix(itotv)==0)disthree(itotv)=1.
                        endif
                    end do
                endif

            endif !if (istep==1) then

            do jiter=1,miter

                rmid=0.
                if(nifsgroup/=0)call load_of_addtional_mass(rmid,deltafi,deltafp)
                call load_of_mass(rmid,deltafi)

                do jstep=1,istep-1
                    beta_mid=dis_shape(:,jstep).d.rmid
                    beta(jstep)=alfa(jstep)*beta_mid
                end do

                rvector=0.0
                do itotv=1,ntotv
                    if(totveq(itotv)/=0.and.totveq(itotv).le.neq_layer1) &
                        rvector(totveq(itotv))=rvector(totveq(itotv))+rmid(itotv)
                end do
                !!int2000
                do itotv=1,ntotv
                    nintf=trans(itotv)%nintf
                    if (nintf/=0) then
                        iieq=totveq(itotv)
                        if(iieq/=0)rvector(iieq)=0.
                        do iintf=1,nintf
                            iieq=totveq(trans(itotv)%listf(iintf))
                            if(iieq/=0) &
                                rvector(iieq)=rvector(iieq)+rmid(itotv)*trans(itotv)%rintf(iintf)
                        end do
                    endif
                end do
                !!int2000


                if (type_solver=='PROFILE') then !zhao 20070829
                    call sksols_layer(global_stiff1,rvector,iseq,1)
                    !             call sksols_layer(global_stiff1,rvector,iseq,0)
                    !!int2000
                    deltafi=0.
                    do itotv=1,ntotv
                        nintf=trans(itotv)%nintf
                        if (iffix(itotv)==0.and.nintf==0) then
                            deltafi(itotv)=rvector(totveq(itotv))
                        elseif(nintf/=0) then
                            listf=>trans(itotv)%listf
                            rintf=>trans(itotv)%rintf
                            do jtotv=1,nintf
                                if(totveq(listf(jtotv))>0)deltafi(itotv)=deltafi(itotv)+rvector(totveq(listf(jtotv)))*rintf(jtotv)
                            enddo
                        endif
                    end do
                    !!int2000
                elseif(type_solver=='SSORPBCG')then !zhao 20070829

                    operation='SOLVE'
                    call solve
                    deltafi=result
                elseif(type_solver=='PARDISO')then
                    operation='SOLVE'
                    call solve
                    deltafi=result
                else
                    stop 'type_solver, in response_spectrum '

                endif
                do jstep=1,istep-1
                    deltafi=deltafi-beta(jstep)*dis_shape(:,jstep)
                end do

                !lamda(istep)=abs(deltafi(max_dofn))
                lamda(istep)=maxval(abs(deltafi))
                !lamda(istep)=0.
                do ipoin=1,0 !npoin
                    do idimn=1,ndimn
                        if(idimn>mdofn)cycle
                        if(lmdofn(idimn)==0)cycle
                        itotv=nodfn(lmdofn(idimn),ipoin)
                        if(itotv==0)cycle
                        if(abs(lamda(istep))<abs(deltafi(itotv)))lamda(istep)=deltafi(itotv)
                    enddo
                enddo
                !maloc=maxloc(abs(deltafi))
                !lamda(istep)=deltafi(maloc(1))

                deltafi=deltafi/lamda(istep)
                print *,'max_dofn=',max_dofn
                print *,'lamda(istep)=',lamda(istep),'lamda_iter=',lamda_iter
                if(abs(lamda(istep)-lamda_iter)/abs(lamda(istep)).le.1.e-12) goto 10
                lamda_iter=lamda(istep)
                print *,'iiter=',jiter,'lamda_iter=',lamda_iter

            end do !for iiter
10          continue
            dis_shape(:,istep)=deltafi

            if (istep==1) then
                loadone=0.
                loadtwo=0.
                loadthree=0.
                if(nifsgroup/=0)call load_of_addtional_mass(loadone,disone,deltafp)
                call load_of_mass(loadone,disone)
                if(nifsgroup/=0)call load_of_addtional_mass(loadtwo,distwo,deltafp)
                call load_of_mass(loadtwo,distwo)
                if(nifsgroup/=0)call load_of_addtional_mass(loadthree,disthree,deltafp)
                call load_of_mass(loadthree,disthree)
            endif

            mdelta=0.
            kdelta=0.
            if(nifsgroup/=0)call load_of_addtional_mass(mdelta,deltafi,deltafp)
            call load_of_mass(mdelta,deltafi)

            alfa_mid=deltafi.d.mdelta
            alfa(istep)=lamda(istep)/alfa_mid

            coefx=deltafi.d.loadone
            coefx=coefx/alfa_mid

            coefy=deltafi.d.loadtwo
            coefy=coefy/alfa_mid

            coefz=deltafi.d.loadthree
            coefz=coefz/alfa_mid

            result_zero=deltafi
            if (nifsgroup/=0)then
                do ipoin=1,npoin
                    idofn=lmdofn(8)
                    itotv=nodfn(idofn,ipoin)
                    if(itotv/=0)result_zero(itotv)=deltafp(itotv)
                end do
            endif
            !only for Mode Superposition Method
            !zhao 2005/06/23
            call load_of_stiff(kdelta,deltafi)

            do itotv=1,ntotv
                write(faiunit)result_zero(itotv)
            enddo

            kstar=deltafi.d.kdelta
            omega=1./sqrt(lamda(istep))
            write(chkunit,100)istep,omega,omega/(2.*3.14159),2.*3.14159/omega,coefx,coefy,coefz,alfa_mid,kstar
100         format('  order=',i5,'  omega=',e12.5,'  freq=',e12.5,'  period=',e12.5,' coefx=',e12.5,' coefy=',e12.5, &
                ' coefz=',e12.5,' mstar=',e12.5,' kstar=',e12.5)

            call residu_f
            if (istep/noutf*noutf==istep)then
                call out_full_write
                if(outplot(1:3)=='GID')   call OUT_GID_WRITE
            endif
            if(istep/nresta*nresta==istep)call resta_read_write(-1)
            if(nforce/=0.or.ngaps/=0)call force_interface
            if((nforce/=0.or.ngaps/=0).and.nextrf==0)call write_force_interface
        end do       !! for istep
    end do     !! loop for iincs

    end subroutine response_spectrum

    SUBROUTINE base_frequency_analysis !20231008

    character(80)text
    integer(ink) itotv,iieq,nintf,iintf,jtotv,  &
        jiter,ipoin,idofn
    integer(ink),pointer::listf(:)
    real   (irk), allocatable::rmid(:),deltafp(:)
    real   (irk),pointer::rintf(:)
    real   (irk) time,lamda_iter,omega,lamda


    call dateandtime(curtime)
    write(chkunit,'(a,a)')'Begin calculating base frequency   at ',curtime
    write(*,'(a,a)')      'Begin calculating base frequency   at ',curtime

    if(allocated(rmid))deallocate(rmid)
    if(allocated(rvector))deallocate(rvector)
    allocate(rmid(ntotv),rvector(ntotv))
    if(nifsgroup/=0)allocate(deltafp(ntotv))

    !    call dfact_time_curve(ttime)
    call modf_var_prescribed

    iiter=1
    call stiff_u
    call mcmatrx('U')
    call hmatrx('W')
    global_stiff1=0.
    call estif_assem_response

    operation='FACTORIZE'
    call solve

    deltafi=0.
    do itotv=1,ntotv
        if(iffix(itotv)==0)then
            deltafi(itotv)=1.
        endif
    enddo

    do jiter=1,100
        rmid=0.
        if(nifsgroup/=0)call load_of_addtional_mass(rmid,deltafi,deltafp)     !nzw 2013-5-6 cancle comment
        call load_of_mass_base(rmid,deltafi)
        rvector=0.0
        do itotv=1,ntotv
            if(totveq(itotv)/=0) then
                rvector(totveq(itotv))=rvector(totveq(itotv))+rmid(itotv)
            endif
        end do

        do itotv=1,ntotv
            nintf=trans(itotv)%nintf
            if(nintf/=0) then
                iieq=totveq(itotv)
                if(iieq/=0)rvector(iieq)=0.
                do iintf=1,nintf
                    iieq=totveq(trans(itotv)%listf(iintf))
                    if(iieq/=0) &
                        rvector(iieq)=rvector(iieq)+rmid(itotv)*trans(itotv)%rintf(iintf)
                end do
            endif
        end do
        operation='SOLVE'
        call solve
        deltafi=result
        lamda=maxval(abs(deltafi))
        deltafi=deltafi/lamda
        if(abs(lamda-lamda_iter)/abs(lamda).le.1.e-8) goto 10
        lamda_iter=lamda

    end do !for iiter
10  continue

    omega=1./sqrt(lamda)
    base_freq=omega
    write(chkunit,100)omega,2.*3.14159/omega
    write(omgunit,100)omega,2.*3.14159/omega     !psy.2018.11.01
100 format('  omega=',e12.5,'  period=',e12.5)

    END SUBROUTINE base_frequency_analysis  !903

    subroutine load_of_mass_base(rmid1,dis) !20231008
    character(1)field1
    integer(ink) ielem,nevab,ic,ielgroup,ii,ipoin,inode,idimn,itotv
    integer(ink),pointer::ldofs(:)
    real   (irk),allocatable::eldis(:),eload(:)
    real   (irk),pointer::fstif(:,:)
    real   (irk) dis(:)
    real(irk) rmid1(:)

    DO igroup =1,ngroup
        field1= group(igroup)%fieldid(1:1)
        if(appear(igroup)>0.and.field1=='U') then
            DO ielgroup = 1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                if(associated(element(ielem)%field(1)%khandmc(2)%fstif)) then
                    fstif=>element(ielem)%field(1)%khandmc(2)%fstif
                    ldofs=>element(ielem)%field(1)%ldofs_f
                    nevab=size(ldofs)
                    allocate(eldis(nevab),eload(nevab))
                    ic=size(fstif,dim=2)
                    eldis = dis(ldofs)
                    if(ic==1) then
                        do ii=1,nevab
                            eload(ii)=fstif(ii,1)*eldis(ii)
                        end do
                    else
                        eload=fstif.x.eldis
                    endif
                    rmid1(ldofs)=rmid1(ldofs)+eload
                end if
                deallocate(eldis,eload)
            end do
        endif
    end do
    !write(7,*)'icaddmass=',icaddmass
    if(icaddmass/=0)then
        do ipoin=1,npoin
            if(icmp(ipoin)==0)cycle
            do idimn=1,ndimn
                itotv=nodfn(idimn,ipoin)
                if(itotv==0)cycle
                rmid1(itotv)=rmid1(itotv)+addmp(idimn,ipoin)*dis(itotv)
            enddo
        enddo
    endif
    end subroutine load_of_mass_base

    subroutine load_of_mass(rmid1,dis)

    character(1)field1
    integer(ink) ielem,nevab,ic,ielgroup,ii,ipoin,inode,idimn,itotv
    integer(ink),pointer::ldofs(:)
    real   (irk),allocatable::eldis(:),eload(:)
    real   (irk),pointer::fstif(:,:)
    real   (irk) dis(:)
    real(irk) rmid1(:)
    DO igroup =1,ngroup
        field1= group(igroup)%fieldid(1:1)
        if (appear(igroup)>0.and.field1=='U') then
            DO ielgroup = 1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                if (associated(element(ielem)%field(1)%khandmc(2)%fstif)) then
                    fstif=>element(ielem)%field(1)%khandmc(2)%fstif
                    ldofs=>element(ielem)%field(1)%ldofs_f
                    nevab=size(ldofs)
                    allocate(eldis(nevab),eload(nevab))
                    ic=size(fstif,dim=2)
                    eldis = dis(ldofs)
                    if (ic==1) then
                        do ii=1,nevab
                            eload(ii)=fstif(ii,1)*eldis(ii)
                        end do
                    else
                        eload=fstif.x.eldis
                    endif
                    rmid1(ldofs)=rmid1(ldofs)+eload
                end if
                deallocate(eldis,eload)
            end do
        endif
    end do

    do ipoin=1,nmcon
        inode=lmcon(ipoin)
        do idimn=1,ndimn
            itotv=nodfn(idimn,inode)
            if (itotv/=0)then
                rmid1(itotv)=rmid1(itotv)+rmcon(idimn,ipoin)*dis(itotv)
                !             write(chkunit,*)ipoin,inode,itotv,rmcon(ipoin),dis(itotv)
            endif
        end do
    end do
    if (icaddmass/=0)then
        do ipoin=1,npoin
            if(icmp(ipoin)==0)cycle
            do idimn=1,ndimn
                itotv=nodfn(idimn,ipoin)
                if(itotv==0)cycle
                rmid1(itotv)=rmid1(itotv)+addmp(idimn,ipoin)*dis(itotv)
            enddo
        enddo
    endif
    end subroutine load_of_mass

    subroutine load_of_stiff(rmid1,dis) !only for Mode Superposition Method

    character(1)field1                  !zhao 2005/06/23
    integer(ink) ielem,nevab,ielgroup,ii,ipoin,inode,idimn,itotv
    integer(ink),pointer::ldofs(:)
    real   (irk),allocatable::eldis(:),eload(:)
    real   (irk),pointer::fstif(:,:)
    real   (irk) dis(:)
    real(irk) rmid1(:)
    DO igroup =1,ngroup
        field1= group(igroup)%fieldid(1:1)
        if (appear(igroup)>0.and.field1=='U') then
            DO ielgroup = 1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                if (associated(element(ielem)%field(1)%khandmc(1)%fstif)) then
                    fstif=>element(ielem)%field(1)%khandmc(1)%fstif
                    ldofs=>element(ielem)%field(1)%ldofs_f
                    nevab=size(ldofs)
                    allocate(eldis(nevab),eload(nevab))
                    eldis = dis(ldofs)
                    eload=fstif.x.eldis
                    rmid1(ldofs)=rmid1(ldofs)+eload
                end if
                deallocate(eldis,eload)
            end do
        endif
    end do

    end subroutine load_of_stiff

    subroutine load_of_addtional_mass(loadmid,dis,deltafp)
    integer(ink) ielem,aelemf,aelems,igroup,jgroup,  &
        idofn,ipea1,ipea2,nnode,itotv
    integer(ink),pointer::ldofs(:)
    real   (irk),pointer::estif(:,:)
    real   (irk),allocatable::value(:),rvectorp(:),eload(:),rvector_mid(:)
    real   (irk) loadmid(:),dis(:),deltafp(:)

    allocate(rvectorp(neq),rvector_mid(neq))
    rvectorp=0. ; rvector_mid=0.

    do ielem=1,nifsgroup
        aelemf=tifs(ielem)%aelemf
        aelems=tifs(ielem)%aelems
        ipea1=0
        igroup=element(aelemf)%group
        if(appear(igroup)>0)ipea1=1
        ipea2=1
        if (aelems/=0) then
            jgroup=element(aelems)%group
            if(appear(jgroup)<=0)ipea2=0
        endif
        if (ipea1==1.and.ipea2==1) then
            ldofs=>tifs(ielem)%ldofs
            estif=>tifs(ielem)%estif
            nnode=size(tifs(ielem)%lnods)
            allocate(value(nnode*ndimn),eload(nnode))
            do idofn=1,nnode*ndimn
                value(idofn)=dis(ldofs(idofn))
            end do
            eload=estif(nnode*ndimn+1:nnode*(ndimn+1),1:nnode*ndimn).x.value
            loadmid(ldofs(nnode*ndimn+1:nnode*(ndimn+1)))=   &
                loadmid(ldofs(nnode*ndimn+1:nnode*(ndimn+1)))+eload
            deallocate(value,eload)
            nullify(ldofs,estif)
        endif
    end do

    rvectorp=0.0
    do itotv=1,ntotv
        if(totveq(itotv)/=0.and.totveq(itotv).gt.neq_layer1) &
            rvectorp(totveq(itotv))=rvectorp(totveq(itotv))+loadmid(itotv)
    end do

    if (type_solver=='PROFILE')THEN
        call sksols_layer(global_stiff1,rvectorp,iseq,2)
    ELSEIF(type_solver=='PARDISO')THEN
        rvector_mid=rvector
        rvector=rvectorp
        operation='SOLVE'
        call solve
        rvectorp=rvector
        rvector=rvector_mid
        !STOP 'type_solver=PROFILE'
    ENDIF
    deltafp=0.
    do itotv=1,ntotv
        if(iffix(itotv)==0.and.totveq(itotv)/=0) &
            deltafp(itotv)=rvectorp(totveq(itotv))
    end do

    loadmid=0.
    do ielem=1,nifsgroup

        aelemf=tifs(ielem)%aelemf
        aelems=tifs(ielem)%aelems
        ipea1=0
        igroup=element(aelemf)%group
        if(appear(igroup)>0)ipea1=1
        ipea2=1
        if (aelems/=0) then
            jgroup=element(aelems)%group
            if(appear(jgroup)<=0)ipea2=0
        endif
        if (ipea1==1.and.ipea2==1) then
            ldofs=>tifs(ielem)%ldofs
            estif=>tifs(ielem)%estif
            nnode=size(tifs(ielem)%lnods)
            allocate(value(nnode),eload(nnode*ndimn))
            do idofn=1,nnode
                value(idofn)=deltafp(ldofs(nnode*ndimn+idofn))
            end do
            eload=estif(1:nnode*ndimn,nnode*ndimn+1:nnode*(ndimn+1)).x.value
            loadmid(ldofs(1:nnode*ndimn))=loadmid(ldofs(1:nnode*ndimn))+eload
            deallocate(value,eload)
            nullify(ldofs,estif)
        endif
    end do

    deallocate(rvectorp,rvector_mid)
10  format(10e15.3)

    end subroutine load_of_addtional_mass
