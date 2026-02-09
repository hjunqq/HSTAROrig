
    PROGRAM FEM90

    !                 M. PASTOR, TONCHUN LI AND P. MIRA
    !                          May, 1997
    !                      All rights reserved.

    use fem_module

    implicit none

    ! Local variables (only used in main program body)
    integer             ie, ig, nel_sub, npoin_sub, ngroup_sub
    integer,allocatable :: list_nel_sub(:), icpoin_sub(:), listpoin_sub(:), &
                           listgroup_sub(:), list_group_sub(:)
    integer,allocatable :: icgroup(:)
    real(irk),allocatable :: coord_sub(:,:)

    open(inpunit,file='inp')
    read (inpunit,*) text
    read (inpunit,*) restart,relis,sysrelis,ADINA,Uopt_R,gamamax !20231215YL
    read (inpunit,*) text
    read (inpunit,*) probn
    !	adina=0

    !mystatus=0
    !CALL GETARG (1, probn,mystatus(1))
    !CALL GETARG (2, t2,mystatus(2))
    !CALL GETARG (3, t3,mystatus(3))
    !if(mystatus(2)>0)read(t2,*)restart
    !if(mystatus(3)>0)read(t3,*)Runblks
    call TIME(char_time)
    print *, 'time: ', char_time
    write(chkunit,*)'time: ', char_time
    ! analysis process
    call global_data ! set the global data, they will be unchanged in the whole


    if(Uopt_R==1)then  !20210502
        read(vcor_unit,*)text
        read(vcor_unit,*)radiusi,vdirect,centerR(1:2)

        read(vcor_unit,*)text
        read(vcor_unit,*)nvarp_U
        !write(7,*)text,nvarp_U
        allocate(varplist_U(2,nvarp_U))
        do i0=1,nvarp_U
            read(vcor_unit,*)varplist_U(:,i0)
            !write(7,*)i0,varplist_U(:,i0)
        end do
        read(vcor_unit,*)text
        read(vcor_unit,*)nintp_U
        !write(7,*)text,nintp_U
        allocate(intplist_U(3,nvarp_U),rintf_U(2,nvarp_U))
        do i0=1,nintp_U
            read(vcor_unit,*)intplist_U(:,i0)
            !write(7,*)i0,intplist_U(:,i0)
            read(vcor_unit,*)rintf_U(:,i0)
            !write(7,*)rintf_U(:,i0)
        end do
        call modify_coord

    endif   !20210502

    !call vsl_gauss_gen()

    !call MKL_VSL_TEST

    !stop


    len1=len_trim(probn)


    !if (outintr.gt.0) then
    !   open(outint,file=probn(1:len1)//'.oit',RECL=npoin*irk,FORM='BINARY',ACCESS='DIRECT')
    !elseif(outintw.gt.0.and.restart/=0)then
    !   open(outint,file=probn(1:len1)//'.oit',FORM='BINARY',ACCESS='append')
    !elseif(outintw.gt.0.and.restart==0)then
    !   open(outint,file=probn(1:len1)//'.oit',FORM='BINARY')
    !endif

    if(type_problem=='Q')then   !20200220
        if(outintr>0) then
            open(outint,file=probn(1:len1)//'.oit',RECL=npoin*irk,FORM='BINARY',ACCESS='DIRECT')
            !elseif(outintw>0)then   !20230402
            !    open(outint,file=probn(1:len1)//'.oit',RECL=npoin*irk,FORM='BINARY')
        endif
    else
        if(outintw>0.and.restart/=0) then  !20230402
            open(outint,file=probn(1:len1)//'.oit',FORM='BINARY', &
                ACCESS='append')
        elseif(outintw>0.and.restart==0) then
            open(outint,file=probn(1:len1)//'.oit',FORM='BINARY')
        endif
    endif

    if(upliftin/=0) &        !20220409
        open(upliftunit,file=probn(1:len1)//'.upf',FORM='BINARY') !20220625

    print *,'nblks=',nblks
    print *,'Input runblks, =?'
    !read *,runblks
    read (inpunit,*)runblks
    call TIME(char_time)
    print *, 'time: ', char_time
    write(chkunit,*)'time: ', char_time
    ! material set
    call material_set
    ! modify the element libary
    call modf_element_lib  !20221124

    call contact_point_to_point  !!ctt2005

    call link_concrete_and_steel  !20210328
    call link_concrete_and_water_pipe !20210411



    !  write(7,*)'ntotv=',ntotv,'nodfn='
    !do ipoin=1,npoin
    !write(7,1992)ipoin,nodfn(:,ipoin)
    !end do

    !call steel_spring_parameter  !!steel 2006
    ! stiffness for interface of Fluid and solid, absorbing boundary
    call stiff_interface_fluid_solid    !!ifs2000
    call stiff_absorb_fluid             !!ifs2000
    call stiff_absorb_solid             !!ifs2000
    call stiff_ifs2006                  !!ifs2006 zhao, 06/03/29

    allocate(toler_var(mdofn))
    call output_read  !20210803
    iwriten=0
    trstep=0

    allocate(result_zero(ntotv)) ; result_zero=0.
    if(upliftin/=0)allocate(uplift_node(npoin))
    if(outind==-1)  allocate(accq(ndimn,npoin))  !20231113

    if(nbackf/=0)then
        allocate(result_zero_g(ntotv)) ; result_zero_g=0.   !20210706
        allocate(result_zero_e(ntotv)) ; result_zero_e=0.   !20210706
    endif
    if(type_problem/='Q'.and.type_problem/='E')allocate(result_first(ntotv))
    if(type_problem=='F')allocate(result_second(ntotv))
    !result_zero=0.0
    if(allocated(result_first))result_first=0.0
    if(allocated(result_second))result_second=0.0
    if(allocated(torel))               torel=0.0
    if(allocated(toforl))              toforl=0.0
    if(nflow/=0)then
        allocate(flowrate(npoin))
        flowrate=0.
    endif

    if (type_problem/='W')then !freq2006
        if(allocated(tofor))deallocate(tofor,stfor,toforl,toform,delitfi,deltafi) !,deltafi_ssorpbcg)
        if(allocated(torel))deallocate(torel)
        if(allocated(floae))deallocate(fmass,floae,floai,fexta)
        allocate(tofor(ntotv),stfor(ntotv),toforl(ntotv),toform(ntotv))
        tofor=0. ; stfor=0. ; toforl=0. ; toform=0.
        if(ninit/=0.and.kinit==2) allocate(torel(ntotv)) ; torel=0.  !20201121
        allocate(delitfi(ntotv),deltafi(ntotv)) !,deltafi_ssorpbcg(ntotv)) !ssorpbcg
    else
        if(allocated(toforw))deallocate(toforw,stforw)
        allocate(toforw(ntotv),stforw(ntotv))
        if(allocated(tofor))deallocate(delitfi,deltafi)
        !allocate(delitfi(ntotv),deltafi(ntotv))
    endif

    if (any(props(:)%name=='NSTOKS')) then
        allocate(fmass(ntotv),floae(ntotv),floai(ntotv),fexta(ntotv))
    endif


    allocate(ice0(nelem))
    ice0=0

    call gid_output_parameter !only for check


    if (any(props(:)%name=='NSTOKS')) then
        allocate(fmass(ntotv),floae(ntotv),floai(ntotv),fexta(ntotv))
        fmass=0.
        floai=0.
        floae=0.
        fexta=0.
    endif

    allocate(line_load_block(nblks),line_temp_block(nblks))   !!rrr
    line_load_block=0
    line_temp_block=0
    nincs=0
    if (restart==0)then
        lblks=0
        lttime=0.0
        lincs=0
    else
        ninit=0
        call resta_read_write(1)
        !!special for sanxia
        if(outintr>0.and.(lblks==outintr-1))lttime=0.   !20200226
        !!end special for sanxia
        if (restart==2)then
            appear(1:ngroup)=appear_process(1:ngroup,lblks)
            call local_stress
            !if (kstab==0.)then !zhao 2010
            if(nforce/=0.or.ngaps/=0)call force_interface
            if((nforce/=0.or.ngaps/=0).and.nextrf==0)call write_force_interface
            !else
            if(kstab/=0.)call safety_factor
            !end if
            call out_full_write
            if(outplot(1:3)=='GID')   call OUT_GID_WRITE
            if(outplot(1:6)=='COSMOS')call OUT_COSMOS_WRITE
            stop
        endif
        do iblks=1,lblks
            read(mainunit,*)text
            read(mainunit,*)nincs
            do iincs=1,nincs
                read(mainunit,*) i0
                read(mainunit,*) f0
            end do
        end do
    end if !if (restart==0)then

    ttime=lttime
    blks_new=lblks+1
    write(*,*)'lblks=',lblks,' blks_new=',blks_new
    incs_new=1
    write(*,*)'lincs=',lincs,' nincs=',nincs
    if (lincs<nincs) then
        incs_new=lincs+1
        blks_new=lblks
    endif
    lineload=0
    write(*,*)'lblks=',lblks,' blks_new=',blks_new
    if(blks_new/=1)lineload=line_load_block(blks_new-1)
    linet=0
    if(blks_new/=1)linet   =line_temp_block(blks_new-1)

    rewind(mainunit)

    do iblks=1,blks_new-1
        read(mainunit,*)text
        read(mainunit,*)nincs
        do iincs=1,nincs
            read(mainunit,*) i0
            read(mainunit,*) f0
        end do
    end do
    lblks=blks_new-1
    lincs=incs_new-1

    print *,'Bparameter=',Bparameter


    if(Bparameter==-1.or.Bparameter==-2)then  !20231030
        call  parameter_back_analysis_verify_read
        call  parameter_back_analysis_verify

    else if((Bparameter>0.and.Bparameter<=2).and.balgor<=1)then
        call parameter_back_analysis_read
        call observe_back_analysis_read
        call parameter_back_analysis(Npara,Mvalue)

    else if((Bparameter>0.and.Bparameter<=2).and.balgor==2)then
        call trust_region_back_analysis_read
        call observe_back_analysis_read
        call trust_region_back_analysis(Npara,Mvalue)
    else if(Bparameter==3)then

        call matrix_rigid_dis
        call observe_back_analysis_read
        call rigid_dis_back_analysis
        stop
    else if(Bparameter==4)then

        call matrix_nodal_value
        call observe_back_analysis_read

        call nodal_value_back_analysis
        stop
    else
        if(nbackdT==2) &     !20230216
            call observe_back_analysis_read
        call process_analysis
    endif

    if(submodel<0) then  !20230407
        write(chkunit,*)'npoin_L,nelem_L,ngroup_L,nstep'
        write(chkunit,*)'ires_u, ires_rot,ires_v, ires_a, ires_T, ires_Tv, ires_P, ires_Pv, ires_Pa'
        write(chkunit,*)'nelgp'

        write(chkunit,1992) npoin,nelem,ngroup,nstep,miter
        write(chkunit,992) res_u,res_rot,res_v, res_a, res_T, res_Tv, res_P, res_Pv, res_Pa
        write(chkunit,1992)group(1:ngroup)%nelgroup


        if(submodel==-1)then
            print *,'input total groups or elements for submodel analysis: ngroup_sub,nel_sub'
            !输入子模型分析的组数或单元数：ngroup_sub,nel_sub
            read *, ngroup_sub,nel_sub
            if(ngroup_sub==0.and.nel_sub==0) goto 11


            if(ngroup_sub/=0)then
                allocate(list_group_sub(ngroup_sub))
                print *,'input list of groups for submodel analysis list_ngroup_sub='
                read *,list_group_sub

                allocate(listgroup_sub(ngroup))
                listgroup_sub=0
                do ig=1,ngroup_sub
                    igroup=list_group_sub(ig)
                    listgroup_sub(igroup)=ig
                end do

                nel_sub=0
                do ig=1,ngroup_sub
                    igroup=list_group_sub(ig)
                    nel_sub=nel_sub+group(igroup)%nelgroup
                end do
                allocate(list_nel_sub(nel_sub))

                nel_sub=0
                do ig=1,ngroup_sub
                    igroup=list_group_sub(ig)
                    do ie=1,group(igroup)%nelgroup
                        ielem=group(igroup)%list(ie)
                        nel_sub=nel_sub+1
                        list_nel_sub(nel_sub)=ielem
                    end do
                end do

                deallocate(list_group_sub)

            elseif(nel_sub/=0)then

                allocate(list_nel_sub(nel_sub))
                print *,'input list of element numbers for  submodel analysis list_nel_sub='
                read *,list_nel_sub
                allocate(icgroup(ngroup))

                icgroup=0
                do ie=1,nel_sub
                    ielem=list_nel_sub(ie)
                    icgroup(element(ielem)%group)=1
                end do
                ngroup_sub=sum(icgroup)

                allocate(listgroup_sub(ngroup))

                ngroup_sub=0
                do ig=1,ngroup
                    if(icgroup(ig)==0)cycle
                    ngroup_sub=ngroup_sub+1
                    listgroup_sub(ig)=ngroup_sub
                end do

                deallocate(icgroup)
            endif

            write(chkunit,*)'ngroup_sub=',ngroup_sub,'nel_sub=',nel_sub

            allocate(icpoin_sub(npoin),listpoin_sub(npoin))
            icpoin_sub=0
            listpoin_sub=0
            do ie=1,nel_sub
                ielem=list_nel_sub(ie)
                lnods=>element(ielem)%field(1)%lnods_f
                icpoin_sub(lnods)=1
                nullify(lnods)
            end do
            npoin_sub=sum(icpoin_sub)
            allocate(coord_sub(ndimn,npoin_sub))
            npoin_sub=0
            do ipoin=1,npoin
                if(icpoin_sub(ipoin)==0)cycle
                npoin_sub=npoin_sub+1
                coord_sub(:,npoin_sub)=coord(:,ipoin)
                listpoin_sub(ipoin)=npoin_sub
            end do

            write(sub_msh_unit,*)'NODES INFORMATION'
            do ipoin=1,npoin_sub
                write(sub_msh_unit,1991)ipoin,coord_sub(:,ipoin)
            end do
            write(sub_msh_unit,*)'ELEMENTS INFORMATION'
            do ie=1,nel_sub
                ielem=list_nel_sub(ie)
                ig=element(ielem)%group
                igroup=listgroup_sub(ig)
                lnods=>element(ielem)%field(1)%lnods_f
                write(sub_msh_unit,1992)ie,size(lnods),listpoin_sub(lnods),igroup
                nullify(lnods)
            end do

            deallocate(list_nel_sub,coord_sub,icpoin_sub,listpoin_sub,listgroup_sub)

        endif
    endif !20230407

11  call out_record_write   !20230407

    if (rmesh<0)then
        rewind(out_msh)
        write(out_msh,*)'mesh dimension = 3 elemtype quadrilateral nnode = 4'
        write(out_msh,*)'coordinates'
        do ipoin=1,npoin
            write(out_msh,1991)ipoin,coord(:,ipoin)
        end do
        write(out_msh,*)'end coordinates'
        write(out_msh,*)'elements'
        tnegid=0
        do igroup=1,ngroup
            if (appear(igroup)==1)then
                do ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    if (ice0(ielem)==0)then
                        tnegid=tnegid+1
                        write(out_msh,1992)tnegid,element(ielem)%field(1)%lnods_f,igroup
                    endif
                enddo
            endif
        enddo

        if (nelem1>0)then
            do igroup=1,ngroup
                if (appear(igroup)==1)then
                    DO ielgroup = 1,group1(igroup)%nelgroup
                        ielem = group1(igroup)%list(ielgroup)
                        if (jce1(ielem)==0)then
                            tnegid=tnegid+1
                            write(out_msh,1992)tnegid,element1(ielem)%field(1)%lnods_f,igroup+ngroup
                        endif
                    end do
                endif
            enddo
        endif

        if (nelem2>0)then
            do igroup=1,ngroup
                if (appear(igroup)==1)then
                    DO ielgroup = 1,group2(igroup)%nelgroup
                        ielem = group2(igroup)%list(ielgroup)
                        tnegid=tnegid+1
                        write(out_msh,992)tnegid,element2(ielem)%field(1)%lnods_f,igroup+2*ngroup
                    end do
                endif
            enddo
        endif
        write(out_msh,*)'end elements'
    endif

1991 format(i10,3e15.3)
1992 format(i10,10i10)
992 format(20i5)

    if(outplot=='GIDL')then
        call GID_CLOSEPOSTRESULTFILE
    endif

    call TIME(char_time)
    print *, 'time: ', char_time
    write(chkunit,*)'time: ', char_time

    END PROGRAM FEM90
