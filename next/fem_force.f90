    subroutine FORCE_EXTERNAL

    character(10) fieldid,class,name,SPtype
    character(20) type_curve
    integer(ink) iplgroup,nudofn,npload,iedge,itcurve,aelem,nrfields,        &
        ifield,ielem,ipload,ielgroup,ipoin,idofn,jdofn,itotv,ndofn, &
        index,nevab,nnode,type_mass,inode,ic,aelems,ipea1,cdbound,imcon, &
        igaps,igapb,ngroupb,igroupb,ic_inertia,ipairs,npairs   !2017/06
    real   (irk) dfact,preact2,xxxx
    integer(ink),pointer::list(:),ldofe(:),ldofs_f(:),lnods(:)
    real   (irk),pointer::rload(:),tload(:),edload(:),fstif(:,:)
    real   (irk),allocatable::value(:),tt(:),cc(:),loadlocal(:,:,:),forceint(:,:),forcel(:,:)
    real   (irk) dx,ca,ss,t1,t2,t3,timer,coordzi,timer1,timer2
    integer(ink) nextr,i0,i1,iextr,iforce,ipface,idimn,itdis,itveloc,matno
    real   (irk),pointer::cordzfree(:)
    real   (irk),pointer::estif0(:,:),estif(:,:)
    real   (irk),allocatable::dfact1(:),dfact2(:),dfact3(:),dfact4(:) !hxl2006 VIE
    real   (irk),allocatable::value_d(:),value_v(:),eload(:),value_s(:)
    real   (irk),allocatable::sxyz(:),xyz1(:),xyz2(:),dsxyz(:,:),speed(:)
    real   (irk) density,e,nu,alfa,beta,g
    real   (irk)  timer0 !20220105
    integer(ink),pointer::ldofs(:)
    integer(ink),allocatable::ic_inertia_group(:)  !2017/06

    write(7,*)'force_external***,edge_load_group=',edge_load_group
    if(allocated(floae))floae=0. !!nstoks
    tofor=0.0
    if(type_load=='ARCLENGTH')tofor_arclength=0.
    !! add point load increment
    do iplgroup=1,nplgroup
        itcurve=pload(iplgroup)%order_time_curve
        type_curve=tcurves(itcurve)%type_curve


        if (type_curve=='EXTRAPOLATION')then !2004/9/11
            nextr=tcurves(itcurve)%nextr
            allocate(cc(nextr))

            dx=tcurves(itcurve)%dx
            ca=tcurves(itcurve)%ca
            ss=ca*ditime/dx
            t1=(2.-ss)*(1.-ss)*.5
            t2=ss*(2.-ss)
            t3=ss*(ss-1.)*.5
            if (nextr==1)then
                cc(1)=1.
            elseif(nextr==2)then
                cc(1)=2;cc(2)=-1.
            elseif(nextr==5)then
                cc(1)=3.;cc(2)=-3.;cc(3)=1.
            endif
            npload=pload(iplgroup)%npload
            allocate(loadlocal(npload,ndimn,nextr),forceint(ndimn,npoin),forcel(ndimn,2*nextr+1))

            do iextr=1,nextr  !iextr
                allocate(tt(2*iextr+1))
                if (iextr==1)then
                    tt(1)=t1;tt(2)=t2;tt(3)=t3
                elseif(iextr==2)then
                    tt(1)=t1**2;tt(2)=2*t1*t2;tt(3)=2*t1*t3+t2**2
                    tt(4)=2*t2*t3;tt(5)=t3**2
                elseif(iextr==5)then
                    tt(1)=t1**3;tt(2)=3*t2*t1**2;tt(3)=3*t1*t2**2+3*t3*t1**2
                    tt(4)=6*t1*t2*t3+t2**3;tt(5)=3*t1*t3**2+3*t3*t2**2
                    tt(6)=3*t2*t3**2;tt(7)=t3**3
                endif

                forceint=0.
                do iforce=1,nforce
                    npface=surface_force(iforce)%npface
                    do ipface=1,npface
                        ipoin=surface_force(iforce)%list_npface(ipface)
                        forceint(:,ipoin)=surface_force(iforce)%ftfor_ext(:,ipface,iextr)
                    end do
                end do
                do i0=1,npload
                    forcel=0.
                    do i1=1,2*iextr+1
                        forcel(:,i1)=forceint(:,pload(iplgroup)%listep(i0,i1))
                    end do

                    do i1=1,ndimn
                        loadlocal(i0,i1,iextr)=forcel(i1,1:2*iextr+1).d.tt
                    end do
                end do
                deallocate(tt)
            end do   !iextr
            list=>pload(iplgroup)%list
            do i0=1,npload
                ipoin=list(i0)
                do i1=1,ndimn
                    jdofn=lmdofn(i1)
                    if (jdofn/=0) then
                        itotv=nodfn(jdofn,ipoin)
                        tofor(itotv)=tofor(itotv)+(cc.d.loadlocal(i0,i1,:))
                    endif
                end do
            end do
            nullify(list)
            deallocate(loadlocal,forceint,forcel,cc)
        elseif(type_curve=='DISCONTROL')then
            ic=tcurves(itcurve)%ttime_curve(1)
            if(iiter==2)preact1=prescrib(ic)%rdofix
            if (iiter>=3)then
                preact2=prescrib(ic)%rdofix
                if(iiter==3)dfact=(preact1-preact0)**2/(2*(preact1-preact0)-(preact2-preact0))
                if(iiter>3)dfact=(preact1-preact0)*preact4/((preact1-preact0)+preact4-(preact2-preact0))
                preact4=dfact
                dfact=preact0+dfact
                dfact=dfact*tcurves(itcurve)%dfact_curve(1)
            else
                dfact=tcurves(itcurve)%dfact_curve(1)
                dfact=dfact*prescrib(ic)%rdofix
            endif
        else
            dfact=tcurves(itcurve)%dfact
            !write(7,*)'istep=',istep,'dfact=',dfact
        endif

        if (type_curve/='EXTRAPOLATION')then !2004/9/11
            nudofn=pload(iplgroup)%nudofn
            npload=pload(iplgroup)%npload
            list=>pload(iplgroup)%list
            do ipload=1,npload
                ipoin=list(ipload)
                do idofn=1,nudofn
                    jdofn=lmdofn(idofn)
                    if (jdofn/=0) then
                        itotv=nodfn(jdofn,ipoin)
                        tofor(itotv)=tofor(itotv)+(pload(iplgroup)%pxyz(idofn))*dfact
                        if(type_curve=='ARCLENGTH')tofor_arclength(itotv)=tofor_arclength(itotv)+pload(iplgroup)%pxyz(idofn)
                    end if
                end do
            end do
            nullify(list)
        endif    !2004/9/11
    end do
    !! initialize tload
    do ielem=1,nelem
        igroup=element(ielem)%group
        nrfields=group(igroup)%nrfields
        !if (appear(igroup)>0) then  !2017/11/19
        do ifield=1,nrfields
            element(ielem)%field(ifield)%tload=0.0
        end do
        !end if   !2017/11/19
    end do

    if (rmesh>0)then
        do ielem=1,nelem1
            igroup=element1(ielem)%group
            nrfields=group(igroup)%nrfields
            !if (appear(igroup)>0.and.jce1(ielem)/=1) then   !2017/11/19
            do ifield=1,nrfields
                element1(ielem)%field(ifield)%tload=0.0    !2017/11/19
            end do
            !end if
        end do
    endif

    if (rmesh>1)then
        do ielem=1,nelem2
            igroup=element2(ielem)%group
            nrfields=group(igroup)%nrfields
            !if (appear(igroup)>0) then   !2017/11/19
            do ifield=1,nrfields
                element2(ielem)%field(ifield)%tload=0.0
            end do
            !end if     !2017/11/19
        end do
    endif
    !! add beam load
    do ielgroup=1,nbeamload
        edload=>beamload(ielgroup)%edload
        aelem=beamload(ielgroup)%aelem
        itcurve=beamload(ielgroup)%itcurve
        dfact=tcurves(itcurve)%dfact
        element(aelem)%field(1)%tload=element(aelem)%field(1)%tload+edload*dfact
        nullify(edload)
    end do
    !! add plate_water pressure
    !write(7,*)'plateload***='
    do ielgroup=1,nplateload
        edload=>plateload(ielgroup)%edload
        aelem=plateload(ielgroup)%aelem
        itcurve=plateload(ielgroup)%itcurve
        dfact=tcurves(itcurve)%dfact
        element(aelem)%field(1)%tload=element(aelem)%field(1)%tload+edload*dfact
        !write(7,*)aelem,element(aelem)%field(1)%tload

        nullify(edload)
    end do
    !! add edge load increment to tload
    do ielgroup=1,edge_load_group
        iedge=edgeload(ielgroup)%iedge
        itcurve=edgeload(ielgroup)%itcurve

        type_curve=tcurves(itcurve)%type_curve
        !print *,'ielgroup=',ielgroup,'iedge=',iedge,'itvurve=',itcurve,'type_curve=',type_curve
        edload=>edgeload(ielgroup)%edload
        dfact=tcurves(itcurve)%dfact

        if (rmesh==0)then  !!!!!!for rmesh==0
            aelem=edges(iedge)%aelem
            !write(7,*)'iedge=',iedge,'aelem=',aelem,'dfact=',dfact,'edload=',edload

            if(type_curve=='ARCLENGTH')ldofs_f => element(aelem)%field(1)%ldofs_f
            ldofe=>edges(iedge)%ldofe
            ndofn=size(ldofe)

            ! print *,'iedge=',iedge,'aelem=',aelem,'ndofn=',ndofn,'size(tload)=',size(element(aelem)%field(1)%tload)
            !print *,'ldofe=',ldofe
            do idofn=1,ndofn
                jdofn=ldofe(idofn)
                element(aelem)%field(1)%tload(jdofn)=        &
                    element(aelem)%field(1)%tload(jdofn)+edload(idofn)*dfact
                if (type_curve=='ARCLENGTH') then
                    itotv=ldofs_f(jdofn)
                    tofor_arclength(itotv)=tofor_arclength(itotv)+edload(idofn)
                endif
            end do
            if(type_curve=='ARCLENGTH')nullify(ldofs_f)
            nullify(ldofe)
            if (allocated(floae))then  !!nstoks
                igroup=element(aelem)%group
                matno = group(igroup)%matno
                name=props(matno)%name
                if (name=='NSTOKS')then
                    ldofs_f => element(aelem)%field(1)%ldofs_f
                    do idofn=1,ndofn
                        jdofn=ldofe(idofn)
                        itotv=ldofs_f(jdofn)
                        floae(itotv)=floae(itotv)+edload(idofn)*dfact
                    end do
                    nullify(ldofs_f)
                endif
            endif
        else    !for rmesh/=0
            lnods=>edges(iedge)%lnode
            do inode=1,size(lnods)
                do idofn=1,ndimn
                    itotv=nodfn(idofn,lnods(inode))
                    if(itotv>0)tofor(itotv)=tofor(itotv)+edload((inode-1)*ndimn+idofn)*dfact
                end do
            end do
            nullify(lnods)
        endif
        nullify(edload)
    end do
    !! add body force to tload

    write(7,*)'body force'
    do ielem=1,nelem
        igroup=element(ielem)%group
        fieldid=group(igroup)%fieldid
        !if (appear(igroup)>0.and.ice0(ielem)/=1) then   !2017/11/19
        nrfields=element(ielem)%nrfields
        do ifield=1,nrfields
            if (associated(element(ielem)%field(ifield)%rload)) then
                rload=> element(ielem)%field(ifield)%rload
                if (fieldid(ifield:ifield)=='U'.or.fieldid(ifield:ifield)=='W')then
                    itcurve=tcurvegravity(igroup)
                    !print *,'ie=',ielem,'igroup=',igroup,'itcurve=',itcurve
                    if (itcurve/=0) then
                        type_curve=tcurves(itcurve)%type_curve
                        ldofs_f => element(ielem)%field(ifield)%ldofs_f
                        dfact=tcurves(itcurve)%dfact
                        element(ielem)%field(ifield)%tload=                  &
                            element(ielem)%field(ifield)%tload+ rload*dfact
                        !write(7,*)'ie=',ielem,'rload=',rload,'dfact=',dfact
                        if (ifield==1.and.allocated(floae))then  !!nstoks
                            matno = group(igroup)%matno
                            name=props(matno)%name
                            if(name=='NSTOKS') &
                                floae(ldofs_f)=floae(ldofs_f)+rload*dfact
                        endif
                        if(type_curve=='ARCLENGTH')tofor_arclength(ldofs_f)=tofor_arclength(ldofs_f)+rload
                        nullify(ldofs_f)
                    endif
                elseif(fieldid(ifield:ifield)=='T') then
                    element(ielem)%field(ifield)%tload=                  &
                        element(ielem)%field(ifield)%tload+ rload
                endif
                nullify(rload)
            endif
        end do
        !endif   !2017/11/19
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!
    if (rmesh>0)then
        do ielem=1,nelem1
            igroup=element1(ielem)%group
            fieldid=group(igroup)%fieldid
            !if (appear(igroup)>0.and.jce1(ielem)/=1) then  !2017/11/19
            nrfields=element1(ielem)%nrfields
            do ifield=1,nrfields
                if (associated(element1(ielem)%field(ifield)%rload)) then
                    rload=> element1(ielem)%field(ifield)%rload
                    if (fieldid(ifield:ifield)=='U'.or.fieldid(ifield:ifield)=='W')then
                        itcurve=tcurvegravity(igroup)
                        if (itcurve/=0) then
                            type_curve=tcurves(itcurve)%type_curve
                            ldofs_f => element1(ielem)%field(ifield)%ldofs_f
                            dfact=tcurves(itcurve)%dfact
                            element1(ielem)%field(ifield)%tload=                  &
                                element1(ielem)%field(ifield)%tload+ rload*dfact
                            nullify(ldofs_f)
                        endif
                    elseif(fieldid(ifield:ifield)=='T') then
                        element1(ielem)%field(ifield)%tload=                  &
                            element1(ielem)%field(ifield)%tload+ rload
                    endif
                    nullify(rload)
                endif
            end do
            !endif   !2017/11/19
        end do
    endif
    !!!!!!!!!!!!!!!!!!!!!
    if (rmesh>1)then
        do ielem=1,nelem2
            igroup=element2(ielem)%group
            fieldid=group(igroup)%fieldid
            !if (appear(igroup)>0) then  !2017/11/19
            nrfields=element2(ielem)%nrfields
            do ifield=1,nrfields
                if (associated(element2(ielem)%field(ifield)%rload)) then
                    rload=> element2(ielem)%field(ifield)%rload
                    if (fieldid(ifield:ifield)=='U'.or.fieldid(ifield:ifield)=='W')then
                        itcurve=tcurvegravity(igroup)
                        if (itcurve/=0) then
                            type_curve=tcurves(itcurve)%type_curve
                            ldofs_f => element2(ielem)%field(ifield)%ldofs_f
                            dfact=tcurves(itcurve)%dfact
                            element2(ielem)%field(ifield)%tload=                  &
                                element2(ielem)%field(ifield)%tload+ rload*dfact
                            nullify(ldofs_f)
                        endif
                    elseif(fieldid(ifield:ifield)=='T') then
                        element2(ielem)%field(ifield)%tload=element2(ielem)%field(ifield)%tload+ rload
                    endif
                    nullify(rload)
                endif
            end do
            !endif  !2017/11/19
        end do
    endif

    !! internal heat source for unsteady temperature problem

    call assemble_boundt_eload


    !! add inertia force

    write(7,*)'interia force***'
    if (type_problem=='F')then


        ic_inertia=0    !ic_inertia,ic_inertia_group 主要用来识别当接触块体完全张开后，不在施加地震惯性力
        allocate(ic_inertia_group(ngroup))
        ic_inertia_group=1

        do igaps=1,ngaps
            npairs=gaps(igaps)%npairs
            do ipairs=1,npairs
                if(gaps(igaps)%state(ipairs)/=0) ic_inertia=1
            end do
        end do
        if(ic_inertia==0)then
            do igapb=1,ngapb !2017/06
                if(gapb(igapb)%nrdof==0)cycle
                ngroupb=gapb(igapb)%ngroupb
                do igroupb=1,ngroupb
                    ic_inertia_group(gapb(igapb)%listgroupb(igroupb))=0
                end do
            end do
        endif   !2017/06

        DO igroup =1,ngroup
            if(force_process(igroup)==0)cycle !zhao 05/08/05

            write(7,*)'igroup=',igroup,'ic_=',ic_inertia_group(:)
            if(ic_inertia_group(igroup)==0) cycle !2017/06
            !if (appear(igroup)>0) then  !2017/11/19
            ! get information from the group level
            fieldid=group(igroup)%fieldid
            class  =group(igroup)%class
            if (fieldid(1:1)=='U'.and.class=='CO')then
                index    =group(igroup)%index
                ndofn    =group(igroup)%dof(1)%nfdof
                nnode    =elkn(index)%el_field(1)%nnode_f
                nevab    =nnode*ndofn
                allocate(value(nevab))
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    if (associated(element(ielem)%field(1)%khandmc(2)%fstif)) then
                        ! if(istep==nstep) &
                        !write(7,*)'ielem=',ielem,'element(ielem)%field(1)%tload0=',element(ielem)%field(1)%tload


                        fstif=>element(ielem)%field(1)%khandmc(2)%fstif
                        ic=size(fstif,dim=2)
                        value=0.0
                        do inode=1,nnode
                            idofn=(inode-1)*ndofn
                            value(idofn+1:idofn+ndimn)=-fachv
                        end do
                        !write(7,*)'ielem=',ielem,'value=',value,'fstif=',fstif,'value=',value
                        if (ic/=1) then
                            element(ielem)%field(1)%tload=element(ielem)%field(1)%tload+matmul(fstif,value)
                        else
                            do idofn=1,nevab
                                element(ielem)%field(1)%tload(idofn)=                  &
                                    element(ielem)%field(1)%tload(idofn)+fstif(idofn,1)*value(idofn)
                            end do
                        endif
                        !if(istep==nstep) &
                        !write(7,*)'ielem=',ielem,'element(ielem)%field(1)%tload=',element(ielem)%field(1)%tload
                        nullify(fstif)
                    endif        !! for associated
                end do          !! for ielgroup
                deallocate(value)
            end if   !! for 'U' and 'CO'
            !endif       !! for aappear  !2017/11/19
        end do         !! for igroup
        deallocate(ic_inertia_group)   !2017/06

        !hxl2006 VIE
        do ielem=1,nabssgroup
            aelems=tabss(ielem)%aelems
            ipea1=0
            igroup=element(aelems)%group
            if(appear(igroup)>0)ipea1=1
            if (ipea1==1) then
                cdbound=tabss(ielem)%cdbound  !!hxl_l  1,for lateral, 2 for bottom
                lnods=>tabss(ielem)%lnods
                cordzfree=>tabss(ielem)%cordzfree

                !write(7,*)'ie=',ielem,'lnods=',lnods,'cordzfree=',cordzfree
                nnode=size(lnods)
                allocate(dfact1(ndimn),dfact2(ndimn),dfact3(ndimn),dfact4(ndimn))
                allocate(value_d(nnode*ndimn),value_v(nnode*ndimn),value_s(nnode*ndimn))
                allocate(sxyz(ndimn),speed(ndimn),dsxyz(ndimn,ndimn),xyz1(ndimn),xyz2(ndimn))
                value_d=0.;value_v=0.;value_s=0.
                matno=element(tabss(ielem)%aelems)%matno
                SPtype=group(element(tabss(ielem)%aelems)%group)%SPtype
                density=props(matno)%mechanical%solid%density !densxx !

                if(Bparameter/=0.and.props(matno)%mechanical%solid%ie/=0)then !20190810
                    e=xvalue(props(matno)%mechanical%solid%ie)
                else
                    e=props(matno)%mechanical%solid%e !exx !
                endif
                if(Bparameter/=0.and.props(matno)%mechanical%solid%iNu/=0)then
                    Nu=xvalue(props(matno)%mechanical%solid%iNu)
                else
                    Nu=props(matno)%mechanical%solid%Nu !uxx !
                endif

                alfa = e*(1-nu)/((1.+nu)*(1.-2.*nu))
                beta = e*nu/((1.+nu)*(1.-2.*nu))
                if(SPtype=='PS')alfa=e/(1.0-nu**2)
                if(SPtype=='PS')beta=alfa*nu
                G= e/(2.*(1.+nu))
                speed(ndimn)=sqrt(alfa/density) !P波波速
                speed(1:(ndimn-1))=sqrt(g/density) !S波波速
                if (cdbound==2)then           !底边界
                    if(hwdirec==0) then !20220105

                        call dfact_time_curve(ttime)
                        dfact1=0.;dfact2=0.;  sxyz=0.
                        do idimn=1,ndimn
                            itdis=earthquake_curve_d(idimn)     !!hxl_l
                            itveloc=earthquake_curve_v(idimn)   !!hxl_l
                            if(itdis>0) &
                                dfact1(idimn)=tcurves(itdis)%dfact     !入射位移波
                            if(itveloc>0)dfact2(idimn)=2.*tcurves(itveloc)%dfact    !入射速度波 *2？
                        end do
                        do inode=1,nnode
                            do idimn=1,ndimn
                                value_d((inode-1)*ndimn+idimn)=dfact1(idimn)
                                value_v((inode-1)*ndimn+idimn)=dfact2(idimn)
                            end do
                        end do


                    else !20220105
                        !****!20220105
                        do inode=1,nnode
                            timer0=0.  !20220105
                            if(hwdirec>0)then !20220105
                                timer0=(coord(hwdirec,lnods(inode))-hcoord)/speed(hwdirec) !20220105
                            else if(hwdirec<0)then !20220105
                                timer0=(hcoord-coord(hwdirec,lnods(inode)))/speed(-hwdirec) !20220105
                            endif  !20220105
                            timer0=ttime-timer0 !计算时间与输入波传播至当前点的时间之差，即已传播至当前点的时间
                            dfact1=0.;dfact2=0.;  sxyz=0.
                            if (timer0>0.)then
                                call dfact_time_curve(timer0)
                                do idimn=1,ndimn
                                    itdis=earthquake_curve_d(idimn)     !!hxl_l
                                    itveloc=earthquake_curve_v(idimn)   !!hxl_l
                                    if(itdis>0)  dfact1(idimn)=tcurves(itdis)%dfact
                                    if(itveloc>0)dfact2(idimn)=tcurves(itveloc)%dfact
                                end do
                            endif
                            do idimn=1,ndimn
                                value_d((inode-1)*ndimn+idimn)=dfact1(idimn)
                                value_v((inode-1)*ndimn+idimn)=dfact2(idimn)
                            end do
                        end do
                        !****!20220105


                    endif !20220105


                elseif(cdbound==1)then       !侧边界
                    dfact1=0.;dfact2=0.;  sxyz=0.
                    do inode=1,nnode
                        timer0=0.  !20220105
                        if(hwdirec>0)then !20220105
                            timer0=(coord(hwdirec,lnods(inode))-hcoord)/speed(hwdirec) !20220105
                        else if(hwdirec<0)then !20220105
                            timer0=(hcoord-coord(hwdirec,lnods(inode)))/speed(-hwdirec) !20220105
                        endif  !20220105
                        coordzi=coord(ndimn,lnods(inode))
                        do idimn=1,ndimn
                            itdis=earthquake_curve_d(idimn)     !!hxl_l
                            itveloc=earthquake_curve_v(idimn)   !!hxl_l
                            timer1=(coordzi-inpcord)/speed(idimn) !输入波传播至当前点的时间
                            timer1=timer1+timer0 !20220105
                            timer1=ttime-timer1 !计算时间与输入波传播至当前点的时间之差，即已传播至当前点的时间
                            if (timer1>0.)then
                                call dfact_time_curve(timer1)
                                if(itdis>0)  dfact1(idimn)=tcurves(itdis)%dfact
                                if(itveloc>0)dfact2(idimn)=tcurves(itveloc)%dfact
                                if(itveloc>0)sxyz(idimn)=-speed(idimn)*density*tcurves(itveloc)%dfact
                                ! sxyz 入射速度波产生的应力г=-ρ*Cs*V
                            endif
                        end do

                        if (ndimn==2)then ! 由1:nidmn-1个切应力及ndimn法向应力推求入射速度波产生应力张量dsxyz
                            dsxyz(1,1)=beta/alfa*sxyz(2)
                            dsxyz(2,2)=sxyz(2)
                            dsxyz(1,2)=sxyz(1)
                            dsxyz(2,1)=sxyz(1)
                        elseif(ndimn==3)then ! 由1:nidmn-1个切应力及ndimn法向应力推求入射速度波产生应力张量dsxyz
                            dsxyz(1,1)=beta/alfa*sxyz(3)
                            dsxyz(2,2)=dsxyz(1,1)
                            dsxyz(3,3)=sxyz(3)
                            dsxyz(1,2)=0.
                            dsxyz(1,3)=sxyz(1)
                            dsxyz(2,1)=0.
                            dsxyz(2,3)=sxyz(2)
                            dsxyz(3,1)=sxyz(1)
                            dsxyz(3,2)=sxyz(2)
                        endif
                        xyz1=dsxyz.x.tabss(ielem)%rr(ndimn,:) !转换至整体坐标系
                        ! xyz1:整体坐标系下上行波产生的应力

                        dfact3=0.
                        dfact4=0.
                        sxyz=0.
                        do idimn=1,ndimn
                            itdis=earthquake_curve_d(idimn)     !!hxl_l
                            itveloc=earthquake_curve_v(idimn)   !!hxl_l
                            timer2=(cordzfree(inode)-inpcord)/speed(idimn)+(cordzfree(inode)-coordzi)/speed(idimn)
                            timer2=timer2+timer0 !20220105
                            !timer2:入射波传播至顶面的时间+从顶面再传播至当前点的时间
                            timer2=ttime-timer2 !计算时间与入射波从顶面反射至当前点的时间之差
                            if (timer2>0.) then
                                call dfact_time_curve(timer2)
                                if(itdis>0)dfact3(idimn)=tcurves(itdis)%dfact
                                if(itveloc>0)dfact4(idimn)=tcurves(itveloc)%dfact
                                if(itveloc>0)sxyz(idimn)=speed(idimn)*density*tcurves(itveloc)%dfact
                                ! sxyz 入射速度波产生的应力г=ρ*Cs*V
                            endif
                        end do



                        if (ndimn==2)then
                            dsxyz(1,1)=beta/alfa*sxyz(2)
                            dsxyz(2,2)=sxyz(2)
                            dsxyz(1,2)=sxyz(1)
                            dsxyz(2,1)=sxyz(1)
                        elseif(ndimn==3)then
                            dsxyz(1,1)=beta/alfa*sxyz(3)
                            dsxyz(2,2)=dsxyz(1,1)
                            dsxyz(3,3)=sxyz(3)
                            dsxyz(1,2)=0.
                            dsxyz(1,3)=sxyz(1)
                            dsxyz(2,1)=0.
                            dsxyz(2,3)=sxyz(2)
                            dsxyz(3,1)=sxyz(1)
                            dsxyz(3,2)=sxyz(2)
                        endif
                        xyz2=dsxyz.x.tabss(ielem)%rr(ndimn,:)
                        ! xyz1:整体坐标系下下行波产生的应力
                        do idimn=1,ndimn
                            value_d((inode-1)*ndimn+idimn)=(dfact1(idimn)+dfact3(idimn)) !自由场位移波
                            value_v((inode-1)*ndimn+idimn)=(dfact2(idimn)+dfact4(idimn)) !自由场速度波
                            value_s((inode-1)*ndimn+idimn)=xyz1(idimn)+xyz2(idimn)       !自由场应力波
                            !实际上将自由场分为两部分：上行波（即输入波）与下行波（即反射波）
                        end do
                    end do
                endif

                !estif  =Int. (RT NT ρ*Cs N R)
                !estif0 =Int. (RT NT k/(2*rb) N R)
                !eload_s=Int. (NT N)

                !将自由场应力转换为结点荷载累加至总体荷载列阵tofor
                ldofs=>tabss(ielem)%ldofs
                allocate(eload(size(ldofs)))
                estif=>tabss(ielem)%estif
                eload=estif.x.value_v

                tofor(ldofs)=tofor(ldofs)+eload
                estif0=>tabss(ielem)%estif0
                eload=estif0.x.value_d
                tofor(ldofs)=tofor(ldofs)+eload
                nullify(estif)
                estif=>tabss(ielem)%eload_s
                eload=estif.x.value_s
                tofor(ldofs)=tofor(ldofs)+eload
                nullify(ldofs,estif,estif0)
                deallocate(eload)
                deallocate(dfact1,dfact2,dfact3,dfact4,value_d,value_v,value_s)
                deallocate(sxyz,speed,dsxyz,xyz1,xyz2)
            endif
        end do
        !end hxl2006 VIE

    end if   !! for fast problems !if (type_problem=='F')then

    !************************************************************************
    !levelset
    !add dynamic pressure of levelset , for fsi problem

    If (type_problem=='F'.and.level_set_problem==2)then

        do iedge=1,nedgel !iedge
            ic=edgesl(iedge)%ic

            if(ic/=1)cycle !ic=1--FSI boundary, ic/=1--others

            edload=>edgesl(iedge)%edload
            aelem=edgesl(iedge)%selem

            !if(type_curve=='ARCLENGTH')  &
            !ldofs_f => element(aelem)%field(1)%ldofs_f

            !dfact=tcurves(itcurve)%dfact
            dfact=1.0
            ldofe=>edgesl(iedge)%ldofe
            ndofn=size(ldofe)

            do idofn=1,ndofn
                jdofn=ldofe(idofn)
                element(aelem)%field(1)%tload(jdofn)=        &
                    element(aelem)%field(1)%tload(jdofn)+edload(idofn)*dfact
            enddo
        enddo
    Endif

    !************************************************************************

    !! add tload to tofor
    do ielem=1,nelem
        igroup=element(ielem)%group
        !if (appear(igroup)>0.and.ice0(ielem)/=1) then    !!2017/11/19

        nrfields=element(ielem)%nrfields
        do ifield=1,nrfields
            if (associated(element(ielem)%field(ifield)%tload)) then
                tload=> element(ielem)%field(ifield)%tload
                ldofe=>element(ielem)%field(ifield)%ldofs_f
                ndofn=size(ldofe)
                do idofn=1,ndofn
                    tofor(ldofe(idofn))=tofor(ldofe(idofn))+tload(idofn)
                    !write(7,*)'ie=',ielem,'itotv=',ldofe(idofn),'tofor=', tofor(ldofe(idofn)),'tload=',tload(idofn)
                end do
                nullify(tload,ldofe)
            end if
        end do
        !endif               !!2017/11/19
    end do

    !if(istep==nstep)then
    !write(7,*)'tofor11='
    !do itotv=1,ntotv
    !write(7,*)itotv,tofor(itotv)
    !end do
    !endif

    !ifs2006 zhao, 06/03/29, icaddmass
    if (icaddmass/=0)then
        do ipoin=1,npoin
            if(icmp(ipoin)==0)cycle
            do idimn=1,ndimn
                itotv=nodfn(idimn,ipoin)
                if(itotv==0)cycle
                xxxx=fachv(idimn)
                tofor(itotv)=tofor(itotv)-addmp(idimn,ipoin)*xxxx
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
                xxxx=fachv(idimn)
                tofor(itotv)=tofor(itotv)-rmcon(idimn,imcon)*xxxx
            enddo
        enddo
    endif

    if(alfa_p4>0)then  !20221124 对应于局部坐标作未知量的节点，将外载进行转换
        allocate(value(ndimn))
        do ipoin=1,npoin
            if (local_p4(ipoin)==0)cycle
            value=0.
            do idofn=1,ndimn
                itotv=nodfn(idofn,ipoin)
                if (itotv/=0)value(idofn)=tofor(itotv)
            end do

            value=prot(:,:,ipoin).x.value
            do idofn=1,ndimn
                itotv=nodfn(idofn,ipoin)
                if (itotv/=0)tofor(itotv)=value(idofn)
            end do

            value=0.
            do idofn=4,2*ndimn
                itotv=nodfn(lmdofn(idofn),ipoin)
                if (itotv/=0)value(idofn-3)=tofor(itotv)
            end do
            value=prot(:,:,ipoin).x.value
            do idofn=1,ndimn   !4,2*ndimn
                itotv=nodfn(ndimn+idofn,ipoin)
                if (itotv/=0)tofor(itotv)=value(idofn)
            end do
        end do
        deallocate(value)
    end if !20221124



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (rmesh>0.and.nelem1>0)then
        do ielem=1,nelem1
            igroup=element1(ielem)%group
            !if (appear(igroup)>0.and.jce1(ielem)/=1) then    !    !2017/11/19
            nrfields=element1(ielem)%nrfields
            do ifield=1,nrfields
                if (associated(element1(ielem)%field(ifield)%tload)) then
                    tload=> element1(ielem)%field(ifield)%tload
                    ldofe=>element1(ielem)%field(ifield)%ldofs_f
                    ndofn=size(ldofe)
                    do idofn=1,ndofn
                        tofor(ldofe(idofn))=tofor(ldofe(idofn))+tload(idofn)
                    end do
                    nullify(tload,ldofe)
                end if
            end do
            !endif               !!2017/11/19
        end do
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (rmesh>1.and.nelem2>0)then
        do ielem=1,nelem2
            igroup=element2(ielem)%group
            !if (appear(igroup)>0) then    !    !2017/11/19
            nrfields=element2(ielem)%nrfields
            do ifield=1,nrfields
                if (associated(element2(ielem)%field(ifield)%tload)) then
                    tload=>element2(ielem)%field(ifield)%tload
                    ldofe=>element2(ielem)%field(ifield)%ldofs_f
                    ndofn=size(ldofe)
                    do idofn=1,ndofn
                        tofor(ldofe(idofn))=tofor(ldofe(idofn))+tload(idofn)
                    end do
                    nullify(tload,ldofe)
                end if
            end do
            !endif               !!2017/11/19
        end do
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(ground_inf/=0.and.allocated(load_space)) then
        tofor(ldofs_space)=tofor(ldofs_space)+load_space
    endif

1111 format(a10,i10,2(a10,f15.5),a10,2f15.5)
1112 format(a10,i10,a10,f15.5,a10,2f15.5)

    end subroutine FORCE_EXTERNAL

    subroutine FORCE_EXTERNAL_w !freq2006

    integer(ink) nnode,inode,cdbound,idimn
    integer(ink),pointer::list(:),ldofe(:),lnods(:)
    complex   (irk),allocatable::dfact1(:),dfact2(:),estif0(:,:),estif(:,:)
    complex   (irk),allocatable::value_d(:),value_v(:),eload(:)
    real   (irk) density,e,nu,alfa,beta,g
    integer(ink),pointer::ldofs(:)


    toforw=0.0
    do ielem=1,nabssgroup
        cdbound=tabss(ielem)%cdbound  ! 1,for lateral, 2 for bottom
        lnods=>tabss(ielem)%lnods
        nnode=size(lnods)
        allocate(dfact1(ndimn),dfact2(ndimn))
        allocate(value_d(nnode*ndimn),value_v(nnode*ndimn))
        value_d=0.;value_v=0.
        if (cdbound==2)then           !bottom
            dfact1=0.;dfact2=0.
            do idimn=1,ndimn
                if (fachv(idimn)/=0)then
                    dfact1(idimn)=cmplx(1.,0.)    !rusheweiyi
                    dfact2(idimn)=2.*cmplx(0.,-ttime)    !rushesudu
                endif
            end do
            do inode=1,nnode
                do idimn=1,ndimn
                    value_d((inode-1)*ndimn+idimn)=dfact1(idimn)
                    value_v((inode-1)*ndimn+idimn)=dfact2(idimn)
                end do
            end do
            ldofs=>tabss(ielem)%ldofs
            allocate(eload(size(ldofs)),estif(size(ldofs),size(ldofs)),estif0(size(ldofs),size(ldofs)))
            estif=cmplx(1.,0.)*tabss(ielem)%estif
            eload=matmul(estif,value_v)    !add
            toforw(ldofs)=toforw(ldofs)+eload

            estif0=cmplx(1.,0.)*tabss(ielem)%estif0
            eload=matmul(estif0,value_d)
            toforw(ldofs)=toforw(ldofs)+eload
            nullify(ldofs)
            deallocate(eload,estif,estif0)
            deallocate(dfact1,dfact2,value_d,value_v)
        endif
    end do

    end subroutine FORCE_EXTERNAL_w

    subroutine assemble_boundt_estif

    character(10)fieldid
    integer(ink) ielgroup,iedge,aelem,ndofn,igroup,nrfields,idofn,jdofn,kdofn,ldofn,ifield
    integer(ink),pointer::ldofe(:)
    real(irk),pointer::edstif(:,:)
    do ielgroup=1,ntelgroup
        iedge=tedgeload(ielgroup)%iedge
        edstif=>tedgeload(ielgroup)%edstif
        aelem=tedges(iedge)%aelem
        ldofe=>tedges(iedge)%ldofe
        ndofn=size(ldofe)
        igroup=element(aelem)%group
        nrfields=group(igroup)%nrfields
        fieldid=group(igroup)%fieldid
        do ifield=1,nrfields
            if(fieldid(ifield:ifield)=='T')goto 1
        end do
1       do idofn=1,ndofn
            do jdofn=1,ndofn
                kdofn=ldofe(idofn)
                ldofn=ldofe(jdofn)
                element(aelem)%field(ifield)%khandmc(1)%fstif(kdofn,ldofn)=     &
                    element(aelem)%field(ifield)%khandmc(1)%fstif(kdofn,ldofn)+     &
                    edstif(idofn,jdofn)
            end do
        end do
        nullify(ldofe,edstif)
    end do

    end subroutine assemble_boundt_estif

    subroutine assemble_boundt_eload

    character(10)fieldid
    integer(ink) ielgroup,iedge,aelem,ndofn,igroup,nrfields,  &
        idofn,jdofn,ifield,itcurve
    integer(ink),pointer::ldofe(:)
    real(irk),pointer::edload(:)
    real(irk) dfact

    print *,'ntelgroup=',ntelgroup
    do ielgroup=1,ntelgroup
        iedge=tedgeload(ielgroup)%iedge
        edload=>tedgeload(ielgroup)%edload
        !      write(chkunit,*)'ielgroup=',ielgroup,'edload=',edload
        aelem=tedges(iedge)%aelem
        ldofe=>tedges(iedge)%ldofe
        ndofn=size(ldofe)
        igroup=element(aelem)%group
        nrfields=group(igroup)%nrfields
        fieldid=group(igroup)%fieldid
        itcurve=tedgeload(ielgroup)%itcurve
        dfact=tcurves(itcurve)%dfact
        do ifield=1,nrfields
            if(fieldid(ifield:ifield)=='T')exit
        end do
        do idofn=1,ndofn
            jdofn=ldofe(idofn)
            element(aelem)%field(ifield)%tload(jdofn)=     &
                element(aelem)%field(ifield)%tload(jdofn)+     &
                edload(idofn)*dfact
        end do
        nullify(ldofe,edload)
    end do

    end subroutine assemble_boundt_eload

    subroutine eload_ifs2006

    integer(ink) iedge,bkind,felem,igroup,nevab,nevabs,nevabf,nnode,inode,idimn,selem,jgroup,idofn,ndofn
    integer(ink) i0,i1,xdir,zdir,jpoin,npseczx,npsecxz,nsect  !20220330
    integer(ink),pointer::ldofs(:),ldofs_s(:),ldofs_f(:)
    real   (irk),pointer::matrix(:,:)
    real   (irk),allocatable::value(:),values(:),valuef(:)
    real   (irk) coef,accx,accz

    if(Icaddmass==3)then  !计算渡槽槽底中心线上水平加速度引起的竖向压力，以及槽底中心线上竖向加速度引起的侧向压力
        !20220330
        igroup=dwpre_aqu%aqu_group
        if(appear(igroup)==0) goto 10
        xdir=dwpre_aqu%xdir;zdir=dwpre_aqu%zdir;nsect=dwpre_aqu%nsect
        npseczx=dwpre_aqu%npseczx;npsecxz=dwpre_aqu%npsecxz
        dwpre_aqu%eloadzx=0.;dwpre_aqu%eloadxz=0.
        do i0=1,nsect
            jpoin=dwpre_aqu%jnode(i0)
            idofn=nodfn(xdir,jpoin)
            accx=result_second(idofn)+fachv(xdir)
            idofn=nodfn(zdir,jpoin)
            accz=result_second(idofn)+fachv(zdir)
            do i1=1,npseczx
                idofn=(i0-1)*npseczx+i1
                dwpre_aqu%eloadzx(idofn)=accx*dwpre_aqu%pzx(i1,i0)
            end do
            do i1=1,npsecxz
                idofn=(i0-1)*npsecxz+i1
                dwpre_aqu%eloadxz(idofn)=accz*dwpre_aqu%pxz(i1,i0)
            end do
        end do

10      continue
    endif  !20220330

    coef=-beeta2*ditime**2

    IF(IFSNEDGE==0)RETURN
    do iedge=1,ifsnedge
        ifsedges(iedge)%eload=0.
    enddo

    do iedge=1,ifsnedge !iedge
        bkind=ifsedges(iedge)%bkind
        if(bkind==2)cycle
        felem=ifsedges(iedge)%felem
        igroup=element(felem)%group
        if(appear(igroup)==0)cycle
        ldofs =>ifsedges(iedge)%ldofs
        matrix=>ifsedges(iedge)%matrix
        nevab=size(ldofs)
        allocate(value(nevab))
        if(bkind==1)value=result_second(ldofs)
        if(bkind==3.or.bkind==4)value=result_first(ldofs)
        ifsedges(iedge)%eload=matrix.x.value
        deallocate(value)
        nullify(matrix,ldofs)
    enddo

    do iedge=1,ifsnedge !iedge
        bkind=ifsedges(iedge)%bkind
        nnode=ifsedges(iedge)%nnode
        if(bkind/=2)cycle
        felem=ifsedges(iedge)%felem
        igroup=element(felem)%group
        selem=ifsedges(iedge)%selem
        jgroup=element(selem)%group
        if(appear(igroup)==0.or.appear(jgroup)==0)cycle
        ldofs =>ifsedges(iedge)%ldofs
        ldofs_s=>ifsedges(iedge)%ldofs_s
        ldofs_f=>ifsedges(iedge)%ldofs_f
        matrix=>ifsedges(iedge)%matrix
        nevab=size(ldofs)
        nevabs=size(ldofs_s)
        nevabf=size(ldofs_f)
        allocate(values(nevabs),valuef(nevabf))
        values=result_second(ldofs_s)
        !      do inode=1,nnode
        !         do idimn=1,ndimn
        !            idofn=(inode-1)*ndimn+idimn
        !            values(idofn)=values(idofn)+fachv(idimn)
        !         enddo
        !      enddo
        ndofn=group(jgroup)%dof(1)%nfdof
        do inode=1,nnode
            idofn=(inode-1)*ndofn
            values(idofn+1:idofn+ndimn)=values(idofn+1:idofn+ndimn)+fachv
        end do
        valuef=result_zero(ldofs_f)
        ifsedges(iedge)%eload(1:nevabs)=-transpose(matrix).x.valuef/coef
        ifsedges(iedge)%eload(nevabs+1:nevab)=matrix.x.values
        deallocate(values,valuef)
        nullify(matrix,ldofs,ldofs_s,ldofs_f)
    enddo

    end subroutine eload_ifs2006

    subroutine eload_ifs2006_w !ifs2006

    integer(ink) iedge,bkind,felem,igroup,nevab,nevabs,nevabf,nnode,inode,idimn,selem,jgroup,idofn,ndofn
    integer(ink),pointer::ldofs(:),ldofs_s(:),ldofs_f(:)
    complex(irk),allocatable::value(:),values(:),valuef(:),matrix(:,:)
    complex(irk) coef

    do iedge=1,ifsnedge
        ifsedges(iedge)%eload=0.
    enddo

    do iedge=1,ifsnedge !iedge
        bkind=ifsedges(iedge)%bkind
        if(bkind==2)cycle
        felem=ifsedges(iedge)%felem
        igroup=element(felem)%group
        if(appear(igroup)==0)cycle
        if(bkind==1)coef=cmplx(-1.0,0.)
        if(bkind==2)cycle
        if(bkind==3.or.bkind==4)coef=cmplx(0.,-1.0/ttime)
        ldofs =>ifsedges(iedge)%ldofs
        nevab=size(ldofs)
        allocate(value(nevab),matrix(nevab,nevab))
        matrix=coef*ifsedges(iedge)%matrix
        value=resultw(ldofs)
        stforw(ldofs)=stforw(ldofs)+matmul(matrix,value)
        deallocate(value,matrix)
        nullify(ldofs)
    enddo

    do iedge=1,ifsnedge !iedge
        bkind=ifsedges(iedge)%bkind
        nnode=ifsedges(iedge)%nnode
        if(bkind/=2)cycle
        coef=cmplx(-1.0,0.)
        felem=ifsedges(iedge)%felem
        igroup=element(felem)%group
        selem=ifsedges(iedge)%selem
        jgroup=element(selem)%group
        if(appear(igroup)==0.or.appear(jgroup)==0)cycle
        ldofs =>ifsedges(iedge)%ldofs
        ldofs_s=>ifsedges(iedge)%ldofs_s
        ldofs_f=>ifsedges(iedge)%ldofs_f
        nevab=size(ldofs)
        nevabs=size(ldofs_s)
        nevabf=size(ldofs_f)
        allocate(values(nevabs),valuef(nevabf),matrix(nevabf,nevabs))
        matrix=coef*ifsedges(iedge)%matrix
        values=resultw(ldofs_s)
        valuef=resultw(ldofs_f)
        stforw(ldofs_s)=stforw(ldofs_s)+matmul(transpose(matrix),valuef)
        stforw(ldofs_f)=stforw(ldofs_f)+matmul(matrix,values)
        deallocate(values,valuef,matrix)
        nullify(ldofs,ldofs_s,ldofs_f)
    enddo

    end subroutine eload_ifs2006_w

    subroutine eload_interface_fluid_solid
    integer(ink) ielem,aelemf,aelems,igroup,jgroup,nevab,  &
        idofn,ipea1,ipea2,nnode,idimn
    integer(ink),pointer::ldofs(:)
    real   (irk),pointer::estif(:,:)
    real   (irk),allocatable::value(:)
    real   (irk) coef

    coef=theta1*ditime

    do ielem=1,nifsgroup
        tifs(ielem)%eload=0.
        aelemf=tifs(ielem)%aelemf
        aelems=tifs(ielem)%aelems
        ipea1=0
        igroup=element(aelemf)%group
        if(appear(igroup)>0)ipea1=1
        ipea2=1
        if (aelems/=0) then
            jgroup=element(aelems)%group
            if(appear(igroup)<=0)ipea2=0
        endif
        if (ipea1==1.and.ipea2==1) then
            ldofs=>tifs(ielem)%ldofs
            estif=>tifs(ielem)%estif
            nevab=size(ldofs)
            allocate(value(nevab))
            nnode=size(tifs(ielem)%lnods)
            !do idofn=1,nnode*ndimn
            !value(idofn)=coef*result_second(ldofs(idofn))
            !end do
            !zhao
            do inode=1,nnode
                do idimn=1,ndimn
                    idofn=(inode-1)*ndimn+idimn
                    value(idofn)=coef*(result_second(ldofs(idofn))+fachv(idimn))
                enddo
            enddo

            do idofn=nnode*ndimn+1,nevab
                value(idofn)=result_zero(ldofs(idofn))
            end do
            tifs(ielem)%eload=estif.x.value
            deallocate(value)
            nullify(ldofs,estif)
        endif
    end do

    end subroutine eload_interface_fluid_solid

    subroutine eload_interface_fs_w !freq2006
    integer(ink) ielem,aelemf,aelems,igroup,jgroup,nevab,idofn,ipea1,ipea2,nnode
    integer(ink),pointer::ldofs(:)
    complex(irk),allocatable::value(:),estif(:,:)
    complex(irk) coef

    coef=cmplx(1.,0.)

    do ielem=1,nifsgroup
        aelemf=tifs(ielem)%aelemf
        aelems=tifs(ielem)%aelems
        ipea1=0
        igroup=element(aelemf)%group
        if(appear(igroup)>0)ipea1=1
        ipea2=1
        if (aelems/=0) then
            jgroup=element(aelems)%group
            if(appear(igroup)<=0)ipea2=0
        endif
        if(ipea1==1.and.ipea2==1) then
            ldofs=>tifs(ielem)%ldofs
            nevab=size(ldofs)
            allocate(value(nevab),estif(nevab,nevab))
            estif=coef*tifs(ielem)%estif
            do idofn=1,nevab
                value(idofn)=resultw(ldofs(idofn))
            end do
            stforw(ldofs)=stforw(ldofs)+matmul(estif,value)
            deallocate(value,estif)
            nullify(ldofs)
        endif
    end do

    end subroutine eload_interface_fs_w

    subroutine eload_absorb_fluid
    integer(ink) ielem,aelemf,igroup,nevab,ipea1
    integer(ink),pointer::ldofs(:)
    real   (irk),pointer::estif(:,:)
    real   (irk),allocatable::value(:)
    real   (irk) coef

    coef=-theta1*ditime

    do ielem=1,nabsfgroup
        tabsf(ielem)%eload=0.
        aelemf=tabsf(ielem)%aelemf
        ipea1=0
        igroup=element(aelemf)%group
        if(appear(igroup)>0)ipea1=1
        if (ipea1==1) then
            ldofs=>tabsf(ielem)%ldofs
            estif=>tabsf(ielem)%estif
            nevab=size(ldofs)
            allocate(value(nevab))
            value=coef*result_first(ldofs)
            tabsf(ielem)%eload=estif.x.value
            deallocate(value)
            nullify(ldofs,estif)
        endif
    end do

    end subroutine eload_absorb_fluid

    subroutine eload_absorb_fluid_w !freq2006

    integer(ink) ielem,aelemf,igroup,nevab,ipea1
    integer(ink),pointer::ldofs(:)
    complex(irk),allocatable::value(:),estif(:,:)
    complex(irk) coef

    coef=-cmplx(0.,ttime)/(ttime**2)

    do ielem=1,nabsfgroup
        aelemf=tabsf(ielem)%aelemf
        ipea1=0
        igroup=element(aelemf)%group
        if(appear(igroup)>0)ipea1=1
        if (ipea1==1) then
            ldofs=>tabsf(ielem)%ldofs
            nevab=size(ldofs)
            allocate(value(nevab),estif(nevab,nevab))
            estif=coef*tabsf(ielem)%estif
            value=resultw(ldofs)
            stforw(ldofs)=stforw(ldofs)+matmul(estif,value)
            deallocate(value,estif)
            nullify(ldofs)
        endif
    end do

    end subroutine eload_absorb_fluid_w

    subroutine eload_absorb_solid

    integer(ink) ielem,aelems,igroup,ipea1
    integer(ink),pointer::ldofs(:)
    real   (irk),pointer::estif(:,:),estif0(:,:)
    real   (irk),allocatable::value(:),eload(:)

    do ielem=1,nabssgroup
        tabss(ielem)%eload=0.
        aelems=tabss(ielem)%aelems
        ipea1=0
        igroup=element(aelems)%group
        if(appear(igroup)>0)ipea1=1

        if (ipea1==1) then

            ldofs=>tabss(ielem)%ldofs
            allocate(value(size(ldofs)),eload(size(ldofs)))
            estif=>tabss(ielem)%estif
            value=result_first(ldofs)
            tabss(ielem)%eload=estif.x.value
            estif0=>tabss(ielem)%estif0
            value=result_zero(ldofs)
            eload=estif0.x.value
            tabss(ielem)%eload=tabss(ielem)%eload+eload
            nullify(ldofs,estif,estif0)
            deallocate(value,eload)
        endif
    end do
    end subroutine eload_absorb_solid

    subroutine eload_back_spring  !20150925

    integer(ink) ielem,itotv

    real   (irk) stif_spring

    do ielem=1,nbspring
        itotv=bspring(ielem)%listdof
        stif_spring=bspring(ielem)%spring
        bspring(ielem)%eload=stif_spring*result_zero(itotv)

    end do
    end subroutine eload_back_spring !20150925

    subroutine eload_absorb_solid_w !freq2006

    integer(ink) ielem,aelems,igroup,ipea1
    integer(ink),pointer::ldofs(:)
    complex(irk),allocatable::value(:),estif(:,:),estif0(:,:)
    complex(irk) coef,coef0

    coef=cmplx(0.,-ttime)
    coef0=cmplx(1.,0.)          !

    do ielem=1,nabssgroup
        aelems=tabss(ielem)%aelems
        ipea1=0
        igroup=element(aelems)%group
        if(appear(igroup)>0)ipea1=1

        if (ipea1==1) then

            ldofs=>tabss(ielem)%ldofs
            allocate(value(size(ldofs)),estif(size(ldofs),size(ldofs)),estif0(size(ldofs),size(ldofs)))
            estif=coef*tabss(ielem)%estif
            value=resultw(ldofs)
            stforw(ldofs)=stforw(ldofs)+matmul(estif,value)
            estif0=coef0*tabss(ielem)%estif0
            value=resultw(ldofs)
            stforw(ldofs)=stforw(ldofs)+matmul(estif0,value)

            nullify(ldofs)
            deallocate(value,estif,estif0)

        endif
    end do

    end subroutine eload_absorb_solid_w

    subroutine FORCE_INTERNAL

    integer(ink) ifield,ielem,nrfields,itotv,ndofn,iedge, &
        inode,ipoin,np_unode,matno,idimn,felem,imcon,   &       !! stablize
        i0,i1,ipairs,npairs_wc,nline_g_w

    real   (irk) xxxx,coef1
    character(10)fieldid,name                    !! stablize
    integer(ink),pointer::ldofe(:),lnods(:),pairnode_wc(:)  !20210417
    real   (irk),pointer::eload(:),ks(:,:)  !20210417
    real   (irk),allocatable::value(:),heat_wc(:)  !20210417

    stfor=0.0
    if(allocated(floai))floai=0.
    !! add tload to tofor
    do ielem=1,nelem
        igroup=element(ielem)%group
        if (appear(igroup)>0.and.ice0(ielem)/=1) then    !
            fieldid= group(igroup)%fieldid
            nrfields=element(ielem)%nrfields
            do ifield=1,nrfields
                if (associated(element(ielem)%field(ifield)%eload)) then
                    eload=> element(ielem)%field(ifield)%eload
                    ldofe=>element(ielem)%field(ifield)%ldofs_f
                    if (fieldid=='UW'.and.allocated(floai).and.ifield==1)then  !!nstoks
                        matno = group(igroup)%matno
                        name=props(matno)%name
                        if(name=='NSTOKS')floai(ldofe)=floai(ldofe)+eload
                    endif
                    ndofn=size(ldofe)
                    do itotv=1,ndofn
                        stfor(ldofe(itotv))=stfor(ldofe(itotv))+eload(itotv)
                        ! if(ldofe(itotv)==5218) &
                        !write(7,*)'ielem=',ielem,'ifield=',ifield,'idofn=',itotv,'stfor=',stfor(ldofe(itotv)),'eload=',eload(itotv)
                    end do
                    nullify(eload,ldofe)
                end if
            end do
        endif               !
    end do
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
        if(icaddmass==3)then   !20220330

            igroup=dwpre_aqu%aqu_group
            if(appear(igroup)/=0)then
                ldofe=>dwpre_aqu%ldofszx
                eload=>dwpre_aqu%eloadzx
                stfor(ldofe)=stfor(ldofe)+eload
                nullify(ldofe,eload)

                ldofe=>dwpre_aqu%ldofsxz
                eload=>dwpre_aqu%eloadxz
                stfor(ldofe)=stfor(ldofe)+eload
                nullify(ldofe,eload)

            endif
        endif !20220330

    endif
    !2013/4/12

    if (nbspring>0)then  !20150925
        do imcon=1,nbspring
            itotv=bspring(imcon)%listdof
            if(itotv==0)cycle
            stfor(itotv)=stfor(itotv)+bspring(imcon)%eload
        enddo
    endif      !20150925

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


    !!!!1
    if (rmesh>0)then
        do ielem=1,nelem1
            igroup=element1(ielem)%group
            if (appear(igroup)>0.and.jce1(ielem)/=1) then    !

                fieldid= group(igroup)%fieldid
                nrfields=element1(ielem)%nrfields
                do ifield=1,nrfields
                    if (associated(element1(ielem)%field(ifield)%eload)) then
                        eload=>element1(ielem)%field(ifield)%eload
                        ldofe=>element1(ielem)%field(ifield)%ldofs_f
                        !write(7,*)'ie=',ielem,'eload=',eload
                        ndofn=size(ldofe)
                        do itotv=1,ndofn
                            stfor(ldofe(itotv))=stfor(ldofe(itotv))+eload(itotv)
                        end do
                        nullify(eload,ldofe)
                    end if
                end do

            endif               !
        end do
    endif

    !!!!2
    if (rmesh>1)then
        do ielem=1,nelem2
            igroup=element2(ielem)%group
            if (appear(igroup)>0) then    !

                fieldid= group(igroup)%fieldid
                nrfields=element2(ielem)%nrfields
                do ifield=1,nrfields
                    if (associated(element2(ielem)%field(ifield)%eload)) then
                        eload=>element2(ielem)%field(ifield)%eload
                        ldofe=>element2(ielem)%field(ifield)%ldofs_f
                        ndofn=size(ldofe)
                        do itotv=1,ndofn
                            stfor(ldofe(itotv))=stfor(ldofe(itotv))+eload(itotv)
                        end do
                        nullify(eload,ldofe)
                    end if
                end do

            endif               !
        end do
    endif

    !!ifs2000
    do ielem=1,nifsgroup
        ldofe=> tifs(ielem)%ldofs
        eload=> tifs(ielem)%eload
        stfor(ldofe)=stfor(ldofe)+eload
        nullify(ldofe,eload)
    end do
    do ielem=1,nabsfgroup
        ldofe=> tabsf(ielem)%ldofs
        eload=> tabsf(ielem)%eload
        stfor(ldofe)=stfor(ldofe)+eload

        nullify(ldofe,eload)
    end do
    do ielem=1,nabssgroup
        ldofe=> tabss(ielem)%ldofs
        eload=> tabss(ielem)%eload
        stfor(ldofe)=stfor(ldofe)+eload
        !write(7,*)'ielem=',ielem,'eload=',eload,'stfor=',stfor(ldofe)

        nullify(ldofe,eload)
    end do

    !!ifs2000

    !ifs2006 zhao, 06/03/29
    if(icaddmass==0)then    !20231215YL 对附加质量法不需要以下集成
        do iedge=1,ifsnedge
            felem=ifsedges(iedge)%felem
            igroup=element(felem)%group
            if(appear(igroup)==0)cycle
            ldofe=>ifsedges(iedge)%ldofs
            eload=>ifsedges(iedge)%eload
            stfor(ldofe)=stfor(ldofe)+eload
            nullify(ldofe,eload)
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
        stfor(ldofs_space)=stfor(ldofs_space)+eload_space
    endif

    if(nwcpipe/=0)then !20210417
        do i0=1,nwcpipe
            nline_g_w=wc_pipe(i0)%nline_g_w
            coef1= wc_pipe(i0)%iwc

            do i1=1,nline_g_w
                npairs_wc=wc_pipe(i0)%line_g_w(i1)%npairs_wc
                allocate(value(npairs_wc),heat_wc(npairs_wc))
                pairnode_wc=>wc_pipe(i0)%line_g_w(i1)%pairnode_wc
                value=0.;heat_wc=0.
                do ipairs=1,npairs_wc
                    itotv=nodfn(lmdofn(10),pairnode_wc(ipairs))
                    if(itotv/=0) &
                        value(ipairs)=result_zero(itotv)
                end do
                !write(7,*)'value=',value
                Ks=>wc_pipe(i0)%line_g_w(i1)%kmatrix_w
                heat_wc=Ks.x.value
                heat_wc=coef1*heat_wc

                do ipairs=1,npairs_wc
                    itotv=nodfn(lmdofn(10),pairnode_wc(ipairs))
                    if(itotv/=0) &
                        stfor(itotv)=stfor(itotv)+heat_wc(ipairs)
                end do
                deallocate(value,heat_wc)
                nullify(Ks,pairnode_wc)
            end do
        end do

    endif  !20210417


    end subroutine FORCE_INTERNAL

    subroutine force_release

    integer(ink) ifield,ielem,nrfields,itotv,ndofn
    integer(ink),pointer::ldofe(:)
    real   (irk),pointer::eload(:)
    integer(ink),allocatable::id(:)


    !! add tload to torel
    do ielem=1,nelem
        igroup=element(ielem)%group
        if (appear(igroup)==-1.or.appear(igroup)==2) then    !

            !write(chkunit,*)'igroup=',igroup,'ielem=',ielem
            nrfields=element(ielem)%nrfields
            do ifield=1,nrfields
                if (associated(element(ielem)%field(ifield)%eload)) then
                    eload=> element(ielem)%field(ifield)%eload
                    ldofe=>element(ielem)%field(ifield)%ldofs_f
                    ndofn=size(ldofe)
                    do itotv=1,ndofn
                        torel(ldofe(itotv))=torel(ldofe(itotv))+eload(itotv)
                    end do
                    !write(chkunit,*)'eload=',eload
                    nullify(eload,ldofe)
                end if
            end do

        endif               !
    end do

    allocate(id(ntotv))
    id=0

    do ielem=1,nelem
        igroup=element(ielem)%group
        if (appear(igroup)==1) then    !

            nrfields=element(ielem)%nrfields
            do ifield=1,nrfields
                ldofe=>element(ielem)%field(ifield)%ldofs_f
                id(ldofe)=1
                nullify(ldofe)
            end do

        endif               !
    end do

    do itotv=1,ntotv
        if(id(itotv)==0)torel(itotv)=0.
    end do
    deallocate(id)

    end  subroutine FORCE_release

    subroutine ELOAD_INITIALIZE

    integer(ink) ifield,ielem,nrfields,igroup
    character(10) fieldi,fieldid

    do ielem=1,nelem
        igroup=element(ielem)%group
        if (appear(igroup)>0.and.ice0(ielem)/=1) then    !
            if(associated(element(ielem)%rh))element(ielem)%rh=0.

            nrfields=element(ielem)%nrfields
            fieldid=group(igroup)%fieldid
            do ifield=1,nrfields
                fieldi=fieldid(ifield:ifield)
                element(ielem)%field(ifield)%eload=0.0_irk
            end do

        endif               !
    end do

    if (rmesh>0)then
        do ielem=1,nelem1
            igroup=element1(ielem)%group
            if (appear(igroup)>0.and.jce1(ielem)/=1) then    !
                !if(associated(element1(ielem)%rh))element1(ielem)%rh=0.

                nrfields=element1(ielem)%nrfields
                do ifield=1,nrfields
                    element1(ielem)%field(ifield)%eload=0.0_irk
                end do

            endif               !
        end do
    endif

    if (rmesh>1)then
        do ielem=1,nelem2
            igroup=element2(ielem)%group
            if (appear(igroup)>0) then    !
                !if(associated(element2(ielem)%rh))element2(ielem)%rh=0.

                nrfields=element2(ielem)%nrfields
                do ifield=1,nrfields
                    element2(ielem)%field(ifield)%eload=0.0_irk
                end do

            endif               !
        end do
    endif

    end   subroutine ELOAD_INITIALIZE
