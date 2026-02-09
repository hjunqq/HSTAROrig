    subroutine modify_coord  !20210502
    integer i1,i2,i0,ipoin,idk,ipoin1,ipoin2
    real(irk) radius0,xy0(2),xyi(2),f1,f2,alfa0

    if(vdirect/=3)then
        i1=1;i2=2
    else
        i1=2;i2=3
    endif

    do i0=1,nvarp_U
        ipoin=varplist_U(1,i0)
        idk=varplist_U(2,i0)
        if(idk==0)cycle
        if(idk==1)then
            coord(i1,ipoin)=sign(radiusi,coord(i1,ipoin))
        else if(idk==2)then
            xy0(1)=coord(i1,ipoin)-centerR(1)
            xy0(2)=coord(i2,ipoin)-centerR(2)
            radius0=sum(xy0**2)
            radius0=sqrt(radius0)
            alfa0=acos(xy0(1)/radius0)
            xyi(1)=radiusi*cos(alfa0);xyi(2)=-radiusi*sin(alfa0)
            xyi=xyi+centerR
            coord(i1,ipoin)=xyi(1)
            coord(i2,ipoin)=xyi(2)
        endif
    end do

    do i0=1,nintp_U
        ipoin=intplist_U(1,i0)
        ipoin1=intplist_U(2,i0)
        ipoin2=intplist_U(3,i0)
        f1=rintf_U(1,i0)
        f2=rintf_U(2,i0)
        coord(i1,ipoin)=coord(i1,ipoin1)*f1+coord(i1,ipoin2)*f2
        coord(i2,ipoin)=coord(i2,ipoin1)*f1+coord(i2,ipoin2)*f2
        !write(chk_unit,*)'ipoin=',ipoin,'ipoin1,2=',ipoin1,ipoin2,'f1=',f1,'f2=',f2
        !write(chk_unit,*)'coord=',coord(1:2,ipoin)
    end do
    end subroutine modify_coord !20210502

    subroutine update_coord_blarge  !20221102
    integer ipoin,idimn,itotv

    do ipoin=1,npoin
        do idimn=1,ndimn
            itotv=nodfn(idimn,ipoin)
            coord(idimn,ipoin)=coord0(idimn,ipoin)+result_zero(itotv)
        end do
    end do

    end subroutine update_coord_blarge !20221102

    subroutine modify_element_information !20210502

    character(10) nameg
    integer(ink)ielknind,ielem,nnode,ngaus,ikg,ig,id,nr_intrules
    real(irk)   djacb,weigp

    integer(ink),pointer::lnods(:)
    real(irk),pointer::elcod(:,:)

    real(irk),allocatable::shape(:),deriv(:,:),cartd(:,:),posgp(:),xjaci(:,:)

    do ielem=1,nelem
        ielknind=element(ielem)%index
        nr_intrules=elkn(ielknind)%nr_intrules
        igroup=element(ielem)%group

        lnods=>element(ielem)%field(1)%lnods_f
        elcod=>element(ielem)%field(1)%elcod_f


        if  (ielknind/=20.and.ielknind/=21)then  !!new
            do ikg=1,nr_intrules     !!!ikg
                ngaus=elkn(ielknind)%ggaus(ikg)%ngaus
                nnode=elkn(ielknind)%ggaus(ikg)%nnode
                nameg=elkn(ielknind)%ggaus(ikg)%name

                allocate(shape(nnode),deriv(ndimn,nnode),cartd(ndimn,nnode),xjaci(ndimn,ndimn))

                allocate(element(ielem)%egaus(ikg)%gpcod(ndimn,ngaus), &
                    element(ielem)%egaus(ikg)%djacb(ngaus))

                if (nameg(1:4)/='mass')allocate(element(ielem)%egaus(ikg)%cartd(ndimn,nnode,ngaus))
                if (nnode==2)djacb=sqrt(sum((elcod(1:ndimn,2)-elcod(1:ndimn,1))**2))

                !coordinate for gauss points & derivatives
                do ig=1,ngaus !ig

                    shape=elkn(ielknind)%ggaus(ikg)%shape(:,ig)
                    deriv=elkn(ielknind)%ggaus(ikg)%deriv(:,:,ig)
                    weigp=elkn(ielknind)%ggaus(ikg)%weigp(ig)

                    do id=1,ndimn
                        element(ielem)%egaus(ikg)%gpcod(id,ig)=sum (elcod(id,1:nnode)*shape(1:nnode))
                    end do

                    if  (nnode==2)then
                        cartd=deriv/djacb
                    else
                        call jacob(ielem, ndimn, nnode,elcod,deriv,cartd, djacb,xjaci)
                    endif

                    if (nameg/='mass')element(ielem)%egaus(ikg)%cartd(:,:,ig)=cartd
                    element(ielem)%egaus(ikg)%djacb(ig)=djacb*weigp
                end do !ig

                deallocate (shape,deriv,cartd,xjaci)

            end do        !!!ikg
        endif
        nullify(lnods,elcod)
    end do

    end subroutine modify_element_information !20210502

    subroutine modf_element_lib

    character(20)material,field1,special,criteria,name,model
    integer(ink) jgroup,matno,nstre,ngvar,ielgroup,ielem,index,ngaus,ngaus1, &
        order_int, order_int1, order_int2,nevab,aevab,nnode,point2x,point1x, &
        inode,edimn,jndex,ig,idimn,icreep,nr,nrfields,jfield,ifield,icr,point1,point2,ii, &
        kinit_g,uplift_ic,kind_wt,nlinkg  !20221124
    integer(ink),pointer::lnods(:)
    integer(ink),allocatable::ngpoin(:,:)
    real(irk), allocatable::shape(:),cartd(:,:),deriv(:,:),elcod0(:,:),a3(:),a1(:), &
        rxp(:,:),xjaci(:,:),rot(:,:)
    real(irk) djacb,weigp,dis,thickness
    ! define gpvar
    write(7,*)'in modf_element_lib'
    if(type_problem=='WT')then !20220409 用于冰雪冻融
        DO jgroup =1,ngroup
            index = group(jgroup)%index
            ngaus =elkn(index)%ggaus(1)%ngaus
            ngvar=2
            DO ielgroup = 1,group(jgroup)%nelgroup
                ielem = group(jgroup)%list(ielgroup)
                allocate(element(ielem)%field(1)%gpvar(ngvar,ngaus))
            end do
        end do
        return
    endif !20220409

    if(alfa_p4>0)then
        allocate(local_p4(npoin),ngpoin(ngroup,npoin))  !20221124
        local_p4=0 !20221124
        ngpoin=0 !20221124
    endif

    DO jgroup =1,ngroup
        !write(7,*)'jgroup=',jgroup
        field1= group(jgroup)%fieldid
        nrfields= group(jgroup)%nrfields
        special= group(jgroup)%special
        jfield=0
        do ifield=1,nrfields
            if (field1(ifield:ifield)=='T') then
                jfield=ifield
                exit
            endif
        end do

        ! get information from the group level
        index = group(jgroup)%index
        matno = group(jgroup)%matno
        kinit_g=group(jgroup)%kinit_g  !20211214
        write(7,*)'jgroup=',jgroup,'kinit_g=',kinit_g
        name=props(matno)%name
        !uplift_ic=group(jgroup)%uplift_ic  !20220409


        if (field1(1:1)=='U')  then
            thickness  =props(matno)%mechanical%solid%thickness   !2017/04/03
            material=props(matno)%mechanical%solid%material

            !if (index.ne.20.and.index.ne.21) then ! not for beam
            if (index.ne.20.and.index.ne.21.and.index/=25) then ! not for beam !steel 2006

                icreep =props(matno)%mechanical%solid%icreep

                !print *,'matno=',matno,'icreep=',icreep
                if(icreep.ne.0)nr=props(matno)%mechanical%solid%creep%nr
                if(material=='GOODMAN')then
                    group(jgroup)%nstre=ndimn
                    model=props(matno)%mechanical%solid%Goodman%model
                endif
                nstre=group(jgroup)%nstre
                ngvar=nstre+3 ! 20210125
                if (index==22)then
                    nstre=8 !! for plate element
                    ngvar=nstre
                endif

                if (index==26)then  !20230910
                    nstre=3 !! for thin_film element
                    ngvar=nstre
                endif

                ngaus=elkn(index)%ggaus(1)%ngaus

                !write(7,*)'jgroup=',jgroup,'matno=',matno,'material=',material
                if (material=='PLANE_LOWFT')then
                    allocate(a1(ndimn),a3(ndimn),rxp(ndimn,ndimn))
                    point1=props(matno)%mechanical%solid%Plane_lowft%point(1)
                    point2=props(matno)%mechanical%solid%Plane_lowft%point(2)
                    dis=sum((coord(:,point2)-coord(:,point1))**2)
                    dis=sqrt(dis)
                    a1=(coord(:,point2)-coord(:,point1))/dis
                    lnods=>props(matno)%mechanical%solid%Plane_lowft%point
                    call normal_local(lnods,a3)
                    !write(7,*)'a3=',a3
                    call direct_goodman(a3,rxp,ndimn,a1)
                    !write(7,*)'matno=',matno,'rxp(1,:)=',rxp(1,:),'rxp(2,:)=',rxp(2,:)
                    deallocate(a1,a3)
                    nullify(lnods)
                endif



                if (material=='GOODMAN') then
                    if(ndimn==2)ngaus =elkn(1)%ggaus(1)%ngaus
                    if(ndimn==3)ngaus =elkn(5)%ggaus(1)%ngaus
                    allocate(a1(ndimn))
                    a1=0.
                    if(ndimn==3)then
                        point1=props(matno)%mechanical%solid%Goodman%point1
                        point2=props(matno)%mechanical%solid%Goodman%point2
                        if(point1/=0.and.point2/=0)then  !20211031
                            dis=sum((coord(:,point2)-coord(:,point1))**2)
                            dis=sqrt(dis)
                            a1=(coord(:,point2)-coord(:,point1))/dis
                        endif !20211031
                        !write(7,*)'a1=',a1
                    endif

                endif

                nnode =elkn(index)%nnode
                nevab =nnode*ndimn !20231215YL

                if(material=='CLASSICALEP'.or.material=='CONCRETE')ngvar=nstre+3  !nstre+epstn+dlan+yvalue
                if(material=='DUNCANCHANG')ngvar=nstre+5  !nstre+q+s+et+vt+p3 !for judgement of downloading
                if(material=='DUNCANCHANG'.and.icreep==3)ngvar=nstre+7  !nstre+q+s+et+vt+p3+tevf+tetf !for creep
                if(material=='SandPZ'.or.   &  !20220409
                    material=='ClayPZ'.or.material=='SoilPZ') ngvar=2*nstre+1
                !if(material=='SoilPZ') ngvar=2*nstre+1  !20220629
            else if(index.eq.20.or.index.eq.21) then   ! for beam

                ngvar=6*(ndimn-1)
                ngaus=1
                nstre=ngvar

            elseif(index==25)then !steel 2006

                ngvar=2*ndimn
                ngaus=1
                nstre=ngvar
                ngvar=ngvar+5 !ngvar+1--for gaptao, ngvar+2--for gapnorm ,ngvar+3--for Ks, ngvar+4--for steel strain
                !ngvar+5--for state 0-close 1-open , integer it first! for lhg ngvar+6 ic_yty

            endif ! for beam

            group(jgroup)%ngvar=ngvar
            group(jgroup)%nstre=nstre


            ! loop for 1:nelgroup
            DO ielgroup=1,group(jgroup)%nelgroup
                ielem=group(jgroup)%list(ielgroup)
                element(ielem)%field(1)%ngvar_f=ngvar

                !if(index.eq.20.or.index.eq.21)then   !20200116
                ! if(props(matno)%geometry%ipd==1) &
                !element(ielem)%rotation=props(matno)%geometry%rotlg
                ! endif
                if(Blarge/=0.and.(index==20.or.index==21.or.index==22.or.index==26)) then  !20221102
                    allocate(element(ielem)%point_direct(2))
                    element(ielem)%point_direct=0
                endif !20221102

                if(alfa_p4>0.and.(index==22.or.index==26)) then  !20230910
                    lnods=>element(ielem)%field(1)%lnods_f
                    ngpoin(jgroup,lnods)=1
                    nullify(lnods)
                endif !20221124

                allocate(element(ielem)%field(1)%gpvar0(ngvar,ngaus))  !20210125
                allocate(element(ielem)%field(1)%gpvar(ngvar,ngaus),element(ielem)%field(1)%sigz(ngaus))
                allocate(element(ielem)%field(1)%bmatx(nstre,nevab,ngaus)) !20231215YL 存储单元B矩阵

                !if(material=='DUNCANCHANG'.and.uplift_ic/=0)then !20220409
                !if(material=='DUNCANCHANG')then !20220607
                !allocate(element(ielem)%field(1)%gpvar_s(ngvar,ngaus))
                kind_wt=props(matno)%mechanical%solid%kind_wt
                if(kind_wt/=0)then
                    allocate(element(ielem)%field(1)%stran0_s(nstre,ngaus))
                    element(ielem)%field(1)%stran0_s=0.

                    allocate(element(ielem)%field(1)%isatu(ngaus))
                    element(ielem)%field(1)%isatu=0
                endif

                !endif    !20220409


                if(material=='STEEL_SP') allocate(element(ielem)%field(1)%kdiag(3*(ndimn-1))) !20211125


                element(ielem)%field(1)%gpvar0=0.0
                element(ielem)%field(1)%gpvar=0.0
                element(ielem)%field(1)%sigz=0.
                element(ielem)%field(1)%bmatx=0.0 !20231215YL
                if(material=='STEEL_EP')  element(ielem)%field(1)%ep=0. !20211125
                if(material=='STEEL_SP')  element(ielem)%field(1)%kdiag=0. !20211125

                if(material=='DUNCANCHANG'.and.type_problem=='F'.and.gamamax/=0)then !20231215YL
                    allocate(element(ielem)%field(1)%gamamax(ngaus))
                    allocate(element(ielem)%field(1)%gamamax0(ngaus))
                    allocate(element(ielem)%field(1)%gamamax_ini(ngaus))
                    allocate(element(ielem)%field(1)%gamamax_error(ngaus))
                    element(ielem)%field(1)%gamamax=0
                    element(ielem)%field(1)%gamamax0=0
                    element(ielem)%field(1)%gamamax_ini=0
                    element(ielem)%field(1)%gamamax_error=0
                endif !20231215YL


                if (material=='CONCRETE')then
                    icr=props(matno)%mechanical%solid%Concrete%icr
                    if (icr==1)then
                        allocate(element(ielem)%field(1)%rr(ndimn,ndimn,ngaus))
                        element(ielem)%field(1)%rr=0.
                    endif
                    if (icr==2.or.icr==3.or.icr==5.or.icr==6) then !zhao09
                        allocate(element(ielem)%field(1)%strain0(nstre+2,ngaus),element(ielem)%field(1)%strain(nstre+2,ngaus))
                        element(ielem)%field(1)%strain0=0.;element(ielem)%field(1)%strain=0.
                        if(icr==6)element(ielem)%field(1)%gpvar0(nstre+1,:)=1
                        if(icr==6)element(ielem)%field(1)%gpvar (nstre+1,:)=1
                    endif
                endif

                if ((jfield/=0.and.field1/='WT').or.icreep.ne.0) then   !!20220409

                    allocate(element(ielem)%field(1)%stran0(nstre,ngaus))
                    element(ielem)%field(1)%stran0=0.

                    if(icreep==4)then   !20180630
                        allocate(element(ielem)%field(1)%vkstrain0(nstre,ngaus),element(ielem)%field(1)%vkstrain(nstre,ngaus))
                        element(ielem)%field(1)%vkstrain0=0.;element(ielem)%field(1)%vkstrain=0.
                    endif

                end if

                if (jfield/=0.and.field1/='WT') then   !! !20220409
                    allocate(element(ielem)%field(1)%dsig(nstre,ngaus))
                    element(ielem)%field(1)%dsig =0.0
                endif


                if ((jfield/=0.and.field1/='WT').and.icreep.ne.0) then   !!20220409
                    allocate(element(ielem)%field(1)%omega(nstre,ngaus,nr))
                    element(ielem)%field(1)%omega=0.
                endif


                if (kinit_g==2) then
                    allocate(element(ielem)%stres0(nstre,ngaus))
                    element(ielem)%stres0=0.0
                end if



                if (material=='CLASSICALEP') then
                    criteria=props(matno)%mechanical%solid%ClassicalEP%criteria
                    if (criteria=='MCJOINT')then
                        allocate(element(ielem)%rotation(1,ndimn),  &
                            element(ielem)%field(1)%ntstress(2,ngaus))
                        element(ielem)%field(1)%ntstress=0.
                        lnods=>element(ielem)%field(1)%lnods_f
                        call normal_local(lnods,element(ielem)%rotation(1,:))
                        nullify(lnods)
                    endif
                endif
                if (name=='NORMK'.or.name=='NOLINORMK')then
                    allocate(element(ielem)%rotation(1,ndimn))
                    lnods=>element(ielem)%field(1)%lnods_f
                    call normal_local(lnods,element(ielem)%rotation(1,:))
                    nullify(lnods)
                endif

                if (material=='PLANE_LOWFT')then
                    allocate(element(ielem)%rotation(ndimn,ndimn),element(ielem)%field(1)%dmatxd(nstre,nstre,ngaus))   !20130510
                    element(ielem)%field(1)%dmatxd=0.
                    element(ielem)%rotation=rxp
                endif


                if (material/='CLASSICALEP'.and.material/='GOODMAN'.and.name=='CONTACT')then
                    allocate(element(ielem)%rotation(1,ndimn),  &
                        element(ielem)%field(1)%ntstress(2,ngaus))
                    lnods=>element(ielem)%field(1)%lnods_f
                    call normal_local(lnods,element(ielem)%rotation(1,:))
                    nullify(lnods)
                endif
                ! for PZ model 20220409
                if(material=='SandPZ'.or.material=='ClayPZ'.or.material=='SoilPZ')then
                    !if(material=='SoilPZ'.or.material=='SandPZ')then
                    order_int1=elkn(index)%el_field(1)%order_intrules(1)
                    ngaus =elkn(index)%ggaus(order_int1)%ngaus
                    allocate(element(ielem)%egaus(order_int1)%iload(ngaus))
                    allocate(element(ielem)%egaus(order_int1)%iload0(ngaus))
                    allocate(element(ielem)%egaus(order_int1)%vdval(6,ngaus))
                    allocate(element(ielem)%egaus(order_int1)%vdval0(6,ngaus))
                    element(ielem)%egaus(order_int1)%iload=0.0
                    element(ielem)%egaus(order_int1)%iload0=0.0
                    element(ielem)%egaus(order_int1)%vdval=0.0
                    element(ielem)%egaus(order_int1)%vdval0=0.0
                endif
                !end for PZ model	!20220409

                !crack 2006
                if(name=='CRACK')allocate(element(ielem)%field(1)%ntstress(1,ngaus))
                !! contact
                if(index==1.and.material=='ELASTIC_SPRING')               &
                    allocate(element(ielem)%field(1)%state(ngaus))

                if(name=='CONTACT')    &
                    allocate(element(ielem)%field(1)%gapg0(ngaus),         &
                    element(ielem)%field(1)%gapg (ngaus),          &
                    element(ielem)%field(1)%gapn0(nnode),          &
                    element(ielem)%field(1)%gapn (nnode),          &
                    element(ielem)%field(1)%state0(ngaus),         &
                    element(ielem)%field(1)%state(ngaus),          &
                    element(ielem)%field(1)%state1(ngaus),         &  !zhao 05/07/22
                    element(ielem)%field(1)%icftcontact(ngaus),    &
                    element(ielem)%field(1)%natural_thickness(ngaus))
                if(name=='CONTACT')element(ielem)%field(1)%icftcontact=0
                !! end contact

                if (name=='CRACK')then !crack 2006
                    allocate(element(ielem)%field(1)%state(ngaus),element(ielem)%field(1)%state0(ngaus))
                    element(ielem)%field(1)%state='close'
                    element(ielem)%field(1)%state0='close'
                endif
                !! Goodman
                if (material=='GOODMAN') then
                    allocate(element(ielem)%rotation(ndimn,ndimn), &
                        element(ielem)%aera_local(ngaus),a3(ndimn),element(ielem)%evk(ndimn,ngaus))
                    !20231215_YL
                    if(model=='WATERTIGHT')then !20231007 止水
                        jndex=1
                        if(ndimn==3)jndex=5
                        order_int=elkn(jndex)%el_field(1)%order_intrules(1)
                        ngaus=elkn(jndex)%ggaus(order_int)%ngaus

                        allocate(element(ielem)%field(1)%relat_dis_gaus0(ndimn,ngaus),	    	&
                            element(ielem)%field(1)%relat_dis_gaus(ndimn,ngaus),		    &
                            element(ielem)%field(1)%relat_dis_nod0(ndimn,nnode/2),		    &
                            element(ielem)%field(1)%relat_dis_nod(ndimn,nnode/2))
                        element(ielem)%field(1)%relat_dis_gaus0=0.;element(ielem)%field(1)%relat_dis_gaus=0.
                        element(ielem)%field(1)%relat_dis_nod0=0.;element(ielem)%field(1)%relat_dis_nod=0.
                    endif

                    !20231215_YL

                    if(model=='FCM')then
                        element(ielem)%field(1)%gpvar0(nstre+1,:)=1.
                        element(ielem)%field(1)%gpvar(nstre+1,:)=1.
                        do igaus=1,ngaus
                            element(ielem)%evk(:,igaus)=   &
                                props(matno)%mechanical%solid%Goodman%fcmp%kns0
                        end do
                    endif

                    lnods=>element(ielem)%field(1)%lnods_f

                    if(model(1:3)=='FCM') then   !20210125
                        allocate(element(ielem)%field(1)%strain0(ndimn,ngaus),    &
                            element(ielem)%field(1)%strain (ndimn,ngaus))  !20210125
                        element(ielem)%field(1)%strain0=0.
                        element(ielem)%field(1)%strain =0.
                    endif
                    call normal_local(lnods,a3)
                    !call direct(a3,element(ielem)%rotation,ndimn)
                    if(ndimn==3.and.(point1==0.or.point2==0))then  !20211031
                        point2x=lnods(2);point1x=lnods(1)
                        dis=sum((coord(:,point2x)-coord(:,point1x))**2)
                        dis=sqrt(dis)
                        a1=(coord(:,point2x)-coord(:,point1x))/dis
                    endif                !20211031
                    call direct_goodman(a3,element(ielem)%rotation,ndimn,a1)

                    !write(7,*)'jgroup=',jgroup,'ie=',ielem, 'rotation='
                    !do idimn=1,ndimn
                    !write(7,*)element(ielem)%rotation(idimn,:)
                    !end do
                    edimn=ndimn-1
                    !jndex=1  !2017/02/14
                    !if(ndimn==3)jndex=5  !2017/02/14

                    jndex=1
                    if (ndimn==3.and.index==9)jndex=5  !2017/02/14
                    if (ndimn==3.and.index==23)jndex=3  !2017/02/14
                    ngaus1=elkn(jndex)%ggaus(1)%ngaus   !2017/02/14
                    allocate(shape(nnode/2),deriv(edimn,nnode/2),cartd(edimn,nnode/2))
                    allocate(elcod0(edimn,nnode/2),xjaci(edimn,edimn))

                    do inode=1,nnode/2
                        do idimn=1,edimn
                            elcod0(idimn,inode)=element(ielem)%rotation(idimn,:).d.coord(:,lnods(inode))
                        end do
                    end do
                    nullify(lnods)

                    !print *,'ielem=',ielem,'index=',index,'jndex=',jndex,'nnode=',nnode,'ngaus=',ngaus

                    do ig=1,ngaus1
                        !print *,'size(shape)=',size(shape),'size(elkn)=',size(elkn(jndex)%ggaus(1)%shape(:,ig))
                        shape=elkn(jndex)%ggaus(1)%shape(:,ig)
                        deriv=elkn(jndex)%ggaus(1)%deriv(:,:,ig)
                        weigp=elkn(jndex)%ggaus(1)%weigp(ig)
                        call jacob(ielem, edimn, nnode/2,elcod0,deriv,cartd, djacb,xjaci)
                        element(ielem)%aera_local(ig)=djacb*weigp*thickness  !2017/04/03
                    end do
                    !!X,Y,Z ---global axis, x,y,z--local axis
                    !!z is the normal direction of the surface
                    !!    if z/=Y, x=Y*z, y=z*x
                    !!    if z=y,  x=X*z, y=z*x
                    deallocate(deriv,cartd,shape,elcod0,a3,xjaci)
                endif
                !! end of Goodman

                !! for Simo & Rifai element
                if ((index==3.or.index==5.or.index==9.or.index==16.or.index==18).and.special(1:1)=='B')then

                    if (ndimn==2) then
                        nevab=8
                        if(special(2:2)=='A') aevab=2
                        if(special(2:2)=='B') aevab=4
                        if(special(2:2)=='C') aevab=7
                        if(special(2:2)=='D') aevab=11
                        if(index==3)nevab=6
                        if(special(2:2)=='B'.and.index==3) aevab=6
                        if(special(2:2)=='C'.and.index==3) aevab=9
                    else if(ndimn==3) then
                        nevab=24
                        if(special(2:2)=='A') aevab=3
                        if(special(2:2)=='B') aevab=9
                        if(special(2:2)=='C') aevab=24
                        if(special(2:2)=='D') aevab=30
                    endif
                    allocate(element(ielem)%rh(aevab))
                    allocate(element(ielem)%alfa(aevab),element(ielem)%alfa_it(aevab))
                    element(ielem)%alfa=0.
                    element(ielem)%alfa_it=0.
                    if (order_time_mdofn(1)>=1)  then
                        allocate(element(ielem)%alfa_first(aevab))
                        element(ielem)%alfa_first=0.
                    endif
                    if (order_time_mdofn(1)==2)  then
                        allocate(element(ielem)%alfa_second(aevab))
                        element(ielem)%alfa_second=0.
                    endif
                    allocate(element(ielem)%estift(aevab,nevab))
                    allocate(element(ielem)%estifh(aevab,aevab))
                    if(field1(1:2)=='UW')allocate(element(ielem)%qmatxa(aevab,nnode))
                endif
            end do  !ielem
            if(material=='GOODMAN') deallocate(a1)
            if(material=='PLANE_LOWFT')deallocate(rxp)
        endif !for field(1:1)='U'

        !! end for Simo & Rifai element
        !! 11/6/04   ! 单纯渗流场考虑非饱和
        if (field1(1:1)=='W'.and.name(1:6)=='NSSoil') then

            order_int=elkn(index)%el_field(1)%order_intrules(1)
            ngaus=elkn(index)%ggaus(order_int)%ngaus
            DO ielgroup=1,group(jgroup)%nelgroup
                ielem=group(jgroup)%list(ielgroup)
                allocate(element(ielem)%egaus(order_int)%permr(ngaus),  &
                    element(ielem)%egaus(order_int)%csmos(ngaus))
                element(ielem)%egaus(order_int)%permr=1.0
                element(ielem)%egaus(order_int)%csmos=0.
                allocate(element(ielem)%egaus(order_int)%pwatr(ngaus))
                element(ielem)%egaus(order_int)%pwatr=0.
            end do
        endif
        !!
        !!23/2/98
        if (field1(1:2)=='UW'.and.name(1:6)=='NSSoil') then
            order_int=elkn(index)%el_field(2)%order_intrules(1)
            ngaus=elkn(index)%ggaus(order_int)%ngaus
            DO ielgroup=1,group(jgroup)%nelgroup
                ielem=group(jgroup)%list(ielgroup)
                allocate(element(ielem)%egaus(order_int)%permr(ngaus))
                element(ielem)%egaus(order_int)%permr=1.0
                allocate(element(ielem)%egaus(order_int)%pwatr(ngaus),   &
                    element(ielem)%egaus(order_int)%csmos(ngaus),   &
                    element(ielem)%egaus(order_int)%poros(ngaus),   &
                    element(ielem)%egaus(order_int)%satur(ngaus),   &
                    element(ielem)%egaus(order_int)%voide(ngaus))


                order_int1=elkn(index)%el_field(2)%order_intrules(2)  !20220707

                if (order_int1/=order_int)then
                    ngaus =elkn(index)%ggaus(order_int1)%ngaus
                    allocate(element(ielem)%egaus(order_int1)%poros(ngaus),   &
                        element(ielem)%egaus(order_int1)%voide(ngaus),   &
                        element(ielem)%egaus(order_int1)%permr(ngaus),    &
                        element(ielem)%egaus(order_int1)%pwatr(ngaus),   &
                        element(ielem)%egaus(order_int1)%csmos(ngaus),   &
                        element(ielem)%egaus(order_int1)%satur(ngaus))
                endif

                order_int1=elkn(index)%el_field(1)%order_intrules(1)  !20220707
                if (order_int1/=order_int)then
                    ngaus =elkn(index)%ggaus(order_int1)%ngaus
                    allocate(element(ielem)%egaus(order_int1)%poros(ngaus),   &
                        element(ielem)%egaus(order_int1)%satur(ngaus))
                endif


                order_int1=elkn(index)%el_field(1)%order_intrules(2)  !20220707
                if (order_int1/=order_int)then
                    ngaus =elkn(index)%ggaus(order_int1)%ngaus
                    allocate(element(ielem)%egaus(order_int1)%poros(ngaus),   &
                        element(ielem)%egaus(order_int1)%satur(ngaus))
                endif

                order_int1=elkn(index)%couple(1)%intrule_couple(1)
                if (order_int1/=order_int)then
                    ngaus =elkn(index)%ggaus(order_int1)%ngaus
                    allocate(element(ielem)%egaus(order_int1)%satur(ngaus))
                endif

                order_int2=elkn(index)%el_field(1)%order_intrules(1)
                ngaus =elkn(index)%ggaus(order_int2)%ngaus
                if(order_int2/=order_int)then
                    allocate(element(ielem)%egaus(order_int2)%pwatr(ngaus), &
                        element(ielem)%egaus(order_int2)%satur(ngaus))
                endif

                if(material=='SandPZ'.or.material=='ClayPZ'.or.material=='SoilPZ')then
                    allocate(element(ielem)%egaus(order_int2)%iload(ngaus))
                    allocate(element(ielem)%egaus(order_int2)%iload0(ngaus))
                    allocate(element(ielem)%egaus(order_int2)%vdval(6,ngaus))
                    allocate(element(ielem)%egaus(order_int2)%vdval0(6,ngaus))
                endif
            end do !ielem
        endif  !!23/2/98
    end do !jgroup

    !20221124
    if(alfa_p4>0)then
        do ipoin=1,npoin
            if(sum(ngpoin(:,ipoin))/=1)cycle   !20231006 ?
            if(ipp4(ipoin)/=1)cycle
            local_p4(ipoin)=1
            itotv=nodfn(1,ipoin)
            nintf=trans(itotv)%nintf
            if(nintf/=0)local_p4(ipoin)=0
        enddo
        deallocate(ngpoin)
    endif
    !20221124


    end  subroutine modf_element_lib

    subroutine normal_local(lnods,rotation)

    integer(ink) lnods(:)
    integer(ink) index,nnode,edimn,order_int,ngaus,ig,inode,nnode1
    integer(ink),allocatable::lnode(:)
    real   (irk) rotation(:),weigp,aa
    real   (irk),allocatable::elcod(:,:),deriv(:,:),s(:,:),a3(:)
    real   (irk),allocatable::shape(:)
    nnode1=size(lnods)
    index=1
    if(ndimn==3.and.nnode1==8)index=5   !2017/02/14
    if(ndimn==3.and.nnode1==6)index=3   !2017/02/14
    nnode=2
    if(ndimn==3.and.nnode1==8)nnode=4   !2017/02/14
    if(ndimn==3.and.nnode1==6)nnode=3   !2017/02/14
    edimn=ndimn-1
    order_int=elkn(index)%el_field(1)%order_intrules(1)
    ngaus=elkn(index)%ggaus(order_int)%ngaus
    allocate(lnode(nnode),elcod(nnode,ndimn))
    lnode=lnods(1:nnode)
    do inode=1,nnode
        elcod(inode,:)=coord(:,lnode(inode))
    end do
    allocate(shape(nnode),deriv(edimn,nnode))
    allocate(s(ndimn,ndimn),a3(ndimn))
    rotation=0.
    do ig=1,ngaus
        shape=elkn(index)%ggaus(order_int)%shape(:,ig)
        deriv=elkn(index)%ggaus(order_int)%deriv(:,:,ig)
        weigp=elkn(index)%ggaus(order_int)%weigp(ig)
        s(1:edimn,:)=MATMUL(deriv,elcod)
        if((edimn+1).eq.3) then
            s(3,1)=s(1,2)*s(2,3)-s(2,2)*s(1,3)
            s(3,2)=s(1,3)*s(2,1)-s(1,1)*s(2,3)
            s(3,3)=s(1,1)*s(2,2)-s(1,2)*s(2,1)
        else
            s(2,1)=-s(1,2)
            s(2,2)=s(1,1)
        endif
        a3=s(edimn+1,:)**2
        aa=sqrt(sum(a3))
        s(edimn+1,:)=s(edimn+1,:)/aa
        a3=s(edimn+1,:)
        rotation=rotation+a3
    end do
    rotation=rotation/ngaus
    rotation=rotation/sqrt(sum(rotation**2)) !p42010
    deallocate(lnode,elcod,shape,deriv,s,a3)

    end subroutine normal_local

    SUBROUTINE modf_var_prescribed

    character(80) text,type_curve,field1
    integer(ink) idofix,itcurve,ldofix,ipoin,i0,itotv,icdofn,mistep,ifixvar, &
        ic,ifield,nrfields,temp_var_curve,ielgroup,ielem,inode,i1
    integer(ink)  vertical_direction  !20230402

    integer(ink),pointer::ldofs_t(:)
    real   (irk) dfact,f0,time0,temp0,time1,temp1,temp2,dfact1,dfact2,dtemp
    real   (irk), allocatable::midt(:),dmidt(:),value(:)

    integer(ink), pointer::ldofixb(:)    !hxl2006 MIF
    integer(ink) ifixvar0,jfixvar    !20230402
    real   (irk) fixed1,fixed2,bb,val_fix,wpres,xc,corz  !20230402,gamaw

    do idofix=1,ndofix
        itcurve =prescrib(idofix)%itcurve
        dfact   =tcurves(itcurve)%dfact
        ldofix  =prescrib(idofix)%ldofix
        ifixvar =prescrib(idofix)%ifixvar  !20230402
        jfixvar=prescrib(idofix)%jfixvar  !20220304

        type_curve=tcurves(itcurve)%type_curve


        if(type_curve=='EQUINCRE')then !2007/9/28
            fixed(ldofix)=fixed(ldofix)+fincre*prescrib(idofix)%vdofix
            ! write(7,*)'istep=',istep,'fincre=',fincre,'ldofix=',ldofix,'fixed=',fixed(ldofix)
        else if(ifixvar==8.and.jfixvar/=0)then   !20230402
            inode=prescrib(idofix)%nodfix
            if  (jfixvar>0) then
                wpres=(dfact-coord(jfixvar,inode))*gamaw
            else
                wpres=(coord(-jfixvar,inode)-dfact)*gamaw
            end if
            !if(val_fix(ifixnods)<0.)val_fix(ifixnods)=0.
            fixed(ldofix)=wpres
        else if(ifixvar==8.and.jfixvar==0)then  !20230402

            fixed(ldofix)=dfact*gamaw

        else if(ifixvar==10.and.jfixvar/=0)then  !20230402

            !将坝体上下游分为不同区，jfixvar=0时，与通常方法相同，由*.loa中的曲线，根据时间来确定温度值；
            !jfixvar/=0时，jfixvar为指定区域沿不同深度随时间变化曲线，可以分为上游和下游。
            inode=prescrib(idofix)%nodfix
            vertical_direction=temp_surface(jfixvar)%vertical_direction
            corz=coord(vertical_direction,inode)
            call surface_point_temp_find(jfixvar,corz, fixed(ldofix))

        else     !20230402


            fixed(ldofix)=dfact*prescrib(idofix)%vdofix
            !print *,'idofix=',idofix,'fixed=',fixed(ldofix)

        endif




        if (type_problem=='F'.and.ntrans>0)then
            ldofixb=>prescrib(idofix)%ldofixb    !ziyouduzhu  hxl
            bb=1.0/(gamaMIF+1.0)
            if (istep==1)then
                if(ifixvar0_inpb==ifixvar0)then
                    fixed(ldofix)=inpru(ldofix)
                else
                    fixed(ldofix)=inpzi(ldofix)
                end if
            end if
            if (istep==2)then
                if (ifixvar0_inpb==ifixvar0)then
                    fixed(ldofix)=2*bb*(k1.d.disA_1(ldofixb))+inpru(ldofix)
                else
                    fixed(ldofix)=2*bb*(k1.d.disB_1(ldofixb))+inpzi(ldofix)
                end if
            end if
            if (istep>=3)then
                if (ifixvar0_inpb==ifixvar0)then
                    fixed1=k1.d.disA_1(ldofixb)
                    fixed2=k2.d.disA_2(ldofixb)
                    fixed(ldofix)=2*bb*fixed1-(bb**2)*fixed2+inpru(ldofix)
                else
                    fixed1=k1.d.disB_1(ldofixb)
                    fixed2=k2.d.disB_2(ldofixb)
                    fixed(ldofix)=2*bb*fixed1-(bb**2)*fixed2+inpzi(ldofix)
                end if
            end if
        end if
        nullify(ldofixb)
        !      result_zero(ldofix)=fixed(ldofix)

    end do

    print *,'outinp=',outinp
    !if (outinp>0.and.iblks>=outinp) then !20220626
    !   allocate(midt(npoin))
    !  icdofn=lmdofn(8)
    !  !read(outinpunit,*)text
    !  read(outinpunit)midt
    !  do ipoin=1,npoin
    !     itotv=nodfn(icdofn,ipoin)
    !     if (itotv/=0) then
    !        result_zero(itotv)=midt(ipoin)
    !        fixed(itotv)=midt(ipoin)
    !     endif
    !  end do
    !
    !  deallocate(midt)
    !endif

    if(outinp>0.and.iblks>=outinp) then
        allocate(midt(npoin))

        icdofn=lmdofn(8)
        read(outinpunit,*)text
        do i0=1,npoin
            read(outinpunit,*)ipoin,midt(ipoin)
        enddo
        do ipoin=1,npoin
            itotv=nodfn(icdofn,ipoin)
            if(itotv/=0) then
                result_zero(itotv)=midt(ipoin)
                fixed(itotv)=midt(ipoin)
            endif
        end do

        deallocate(midt)
    endif


    if(outintr<0)then  !20200226
        DO igroup =1,ngroup
            if(appear(igroup)>0) then
                field1=group(igroup)%fieldid
                nrfields=group(igroup)%nrfields
                ic=0
                do ifield=1,nrfields
                    if(field1(ifield:ifield)=='T')then
                        ic=1
                        exit
                    end if
                end do
                if(ic==0) goto 1
                time0=group(igroup)%temp_pre%time0
                temp0=group(igroup)%temp_pre%temp0
                temp_var_curve=group(igroup)%temp_pre%temp_var_curve
                group(igroup)%temp_pre%dtemp=0.
                temp1=temp0
                if(ttime>=time0.and.temp_var_curve/=0)then
                    time1=ttime-ditime*inc_step
                    if(time1>=time0)then
                        call dfact_temp_pre(temp_var_curve,time1-time0,dfact1)
                        temp1=temp0+dfact1
                    endif
                    temp2=temp0
                    if(ttime>time0)then
                        call dfact_temp_pre(temp_var_curve,ttime-time0,dfact2)
                        temp2=temp0+dfact2
                    endif
                    group(igroup)%temp_pre%dtemp=temp2-temp1
                    do ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        ldofs_t=>element(ielem)%field(ifield)%ldofs_f
                        deltafi(ldofs_t)=temp2-temp1    !20200226
                        result_zero(ldofs_t)=temp2 !20200226
                        fixed(ldofs_t)=temp2  !20200226
                        nullify(ldofs_t)
                    enddo
                endif
1               continue
            endif  !if(appear(igroup)>0) then
        end do !igroup
    endif !if(outintr<0)then  !20200226

    if(outintr>0.and.iblks>=outintr) then !20200220
        if(.not.allocated(midt))allocate(midt(npoin))
        midt=0.
        icdofn=lmdofn(10)
        if(inc_step==1) then
            read(outint,rec=trstep)midt
        else
            allocate(dmidt(npoin))
            do mistep=trstep-inc_step+1,trstep  !cj042 add step interpolation for thermal stress 20191120
                print *, 'mistep=',mistep,'trstep=',trstep
                read(outint,rec=mistep)dmidt
                midt=midt+dmidt
            end do
            deallocate(dmidt)
        end if
        write(chkunit,*)'temperature**'
        do ipoin=1,npoin
            if (nodfn(icdofn,ipoin).gt.0) then

                deltafi(nodfn(icdofn,ipoin))=midt(ipoin)    !20200226
                result_zero(nodfn(icdofn,ipoin))=result_zero(nodfn(icdofn,ipoin))+midt(ipoin) !20200226
                fixed(nodfn(icdofn,ipoin))=result_zero(nodfn(icdofn,ipoin))  !20200226

                !result_zero(nodfn(icdofn,ipoin))=midt(ipoin)
                !fixed(nodfn(icdofn,ipoin))=midt(ipoin)

            endif
        end do
    endif

    if(upliftin>0.and.iblks>=upliftin)then  !20221119
        if(water_level(iblks)>=0.)then !20230331
            do ipoin=1,npoin
                xc=water_level(iblks)-coord(ndimn,ipoin)
                if(xc<0.)cycle

                uplift_node(ipoin)=xc*gamaw   !20230402
            enddo
        else   !20230331
            !write(7,*)'uplift_node='
            read(upliftunit)uplift_node
            do ipoin=1,npoin
                if(uplift_node(ipoin)<=0.)uplift_node(ipoin)=0.
                !if(uplift_node(ipoin)>0.) &
                !write(7,*)ipoin,uplift_node(ipoin)
            end do
        endif
    endif !20221119

    if(outind==-1) then  !20231113
        accq=0.
        allocate(value(ndimn))
        do i0=1,tbpointsu
            read(outindunit,*)i1,value
            ipoin=listbpointsu_t(i0)
            accq(:,ipoin)=value
        enddo
        deallocate(value)
    endif !20231113


    END SUBROUTINE modf_var_prescribed

    subroutine dfact_temp_pre(itcurve,time,dfact) !20200226

    character(20)type_curve
    real(irk) time,time1,time2,fact1,fact2,dfact,a0sin,asin,wsin,w0sin
    real(irk)theta0,expon,halftime,totime   	! cj042 20191104 equivalent age
    integer(ink) itcurve,itime,ntime

    ntime=tcurves(itcurve)%ntime
    type_curve=tcurves(itcurve)%type_curve

    if(type_curve=='LINEAR'.or.type_curve=='LNLINEAR') then
        if(time<tcurves(itcurve)%ttime_curve(1)) then
            dfact=0.0
        else if(time>=tcurves(itcurve)%ttime_curve(ntime)) then
            dfact=tcurves(itcurve)%dfact_curve(ntime)
        else
            do itime=1,ntime-1
                time1=tcurves(itcurve)%ttime_curve(itime)
                time2=tcurves(itcurve)%ttime_curve(itime+1)
                dfact=0.0
                if(time>=time1.and.time<time2) then
                    fact1=tcurves(itcurve)%dfact_curve(itime)
                    fact2=tcurves(itcurve)%dfact_curve(itime+1)
                    dfact=fact1+(time-time1)/(time2-time1)*(fact2-fact1)
                    exit
                endif
            end do
        end if
        if(type_curve=='LNLINEAR')dfact=exp(dfact)

    elseif(type_curve=='HARMONIC')  then
        a0sin=tcurves(itcurve)%a0sin
        asin =tcurves(itcurve)%asin
        wsin =tcurves(itcurve)%wsin
        w0sin=tcurves(itcurve)%w0sin
        dfact=a0sin+asin*sin(wsin*time+w0sin)

    elseif(type_curve=='DEXPONENTIAL')  then	!d(f)=a*b*exp(b*t)
        time1=tcurves(itcurve)%ttime_curve(1)    !! b
        fact1=tcurves(itcurve)%dfact_curve(1)    !! a
        dfact=fact1*time1*exp(fact1*time)
    elseif(type_curve=='DABT')  then	!d(f)=a*b/(b+t)**2
        time1=tcurves(itcurve)%ttime_curve(1)    !! a
        fact1=tcurves(itcurve)%dfact_curve(1)    !! b
        dfact=fact1*time1/(fact1+time)**2

    else
        print *, 'no type_curve'
        stop
    end if
    end subroutine dfact_temp_pre !20200226

    SUBROUTINE modf_var_prescribed_w !freq2006

    integer(ink) idofix,ldofix
    do idofix=1,ndofix
        ldofix  =prescrib(idofix)%ldofix
        fixed(ldofix)=(prescrib(idofix)%vdofix)/(ttime**2)
    end do

    END SUBROUTINE modf_var_prescribed_w

    subroutine dfact_time_curve(time)

    character(20)type_curve
    real   (irk) time,time1,time2,fact1,fact2,dfact,dtrec,dtend,ample,a0sin,asin,wsin,w0sin
    integer(ink) itcurve,itime,ntime,mgash,ngash,nf,ii,NFs,i
    real   (irk) dtbegin,pai   !hxl
    real (irk) avTem,dTem,detTime,stTime  !20200220
    real   (irk),allocatable::ai(:),omega(:)

    do itcurve=1,ntcurve
        ntime=tcurves(itcurve)%ntime
        type_curve=tcurves(itcurve)%type_curve
        !print *,'itcurve=',itcurve, 'no type_curve=',type_curve
        if (type_curve=='LINEAR'.or.type_curve=='LNLINEAR'.or.type_curve=='WATERLEVEL') then
            if (time<tcurves(itcurve)%ttime_curve(1)) then
                dfact=0.0
            else if(time>=tcurves(itcurve)%ttime_curve(ntime)) then
                dfact=tcurves(itcurve)%dfact_curve(ntime)
            else
                do itime=1,ntime-1
                    time1=tcurves(itcurve)%ttime_curve(itime)
                    time2=tcurves(itcurve)%ttime_curve(itime+1)
                    dfact=0.0
                    if (time>=time1.and.time<time2) then
                        fact1=tcurves(itcurve)%dfact_curve(itime)
                        fact2=tcurves(itcurve)%dfact_curve(itime+1)
                        dfact=fact1+(time-time1)/(time2-time1)*(fact2-fact1)
                        exit
                    endif
                end do
            end if
            if(type_curve=='LNLINEAR')dfact=exp(dfact)

        elseif(type_curve=='TEMPERATURE')  then
            time1=tcurves(itcurve)%ttime_curve(1)
            time2=tcurves(itcurve)%ttime_curve(2)
            fact1=tcurves(itcurve)%dfact_curve(1)
            fact2=tcurves(itcurve)%dfact_curve(2)
            dfact=time1+fact1*sin(2*3.14159*(time2+time)/fact2)
        elseif(type_curve=='COS')  then !17.5+10.8*COS(2*3.14/12*({TIME}/30-6.8)) 20200220
            avTem=tcurves(itcurve)%dfact_curve(1)   !!cj042@126.com  2013.4.18
            dTem=tcurves(itcurve)%dfact_curve(2)
            detTime=tcurves(itcurve)%dfact_curve(3)
            stTime=tcurves(itcurve)%dfact_curve(4)
            dfact=avTem+dTem*cos(2*3.14159/12*((time+stTime)/30.-detTime))

        elseif(type_curve=='PEAK')  then
            dfact=0.0
            do itime=1,ntime
                time1=tcurves(itcurve)%ttime_curve(itime)
                if (time>=time1) then
                    dfact=tcurves(itcurve)%dfact_curve(itime)
                endif
            end do
        elseif(type_curve=='HARMONIC')  then
            a0sin=tcurves(itcurve)%a0sin
            asin =tcurves(itcurve)%asin
            wsin =tcurves(itcurve)%wsin
            w0sin=tcurves(itcurve)%w0sin
            dtbegin=tcurves(itcurve)%dtbegin
            dtend=tcurves(itcurve)%dtend
            if (time<dtbegin)then
                dfact=1.e-30
            else if(time>(dtbegin+dtend)) then
                dfact=1.e-30
            else
                dfact=a0sin+asin*sin(wsin*time+w0sin)
            endif
        elseif(type_curve=='FOURIERSERIES')  then
            Nfs=tcurves(itcurve)%NFS
            !write(7,*)'nfs=',nfs
            dtbegin=tcurves(itcurve)%dtbegin
            dtend=tcurves(itcurve)%dtend
            allocate(ai(nfs),omega(nfs))
            !write(7,*)'ai=',tcurves(itcurve)%ai,'omega=',tcurves(itcurve)%omega
            ai=tcurves(itcurve)%ai
            omega=tcurves(itcurve)%omega

            if (time<dtbegin)then
                dfact=1.e-30
            else if(time>(dtbegin+dtend)) then
                dfact=1.e-30
            else
                pai=3.14159
                dfact=0.
                do i=1,nfs
                    dfact=dfact+ai(i)*sin(omega(i)*pai*time)
                end do
            endif

            write(7,*)'time=',time,'dfact=',dfact
            deallocate(ai,omega)

        elseif (type_curve=='SEISMIC')  then
            dtrec=tcurves(itcurve)%dtrec   !hxl
            dtbegin=tcurves(itcurve)%dtbegin
            dtend=tcurves(itcurve)%dtend
            ample=tcurves(itcurve)%ample
            dfact=0.0
            if (time<dtbegin)then
                dfact=1.e-30
            else if(time>(dtbegin+dtend)) then
                dfact=1.e-30
            else
                time1=(TIME-dtbegin)/DTREC+1
                MGASH=time1
                NGASH=MGASH+1
                time1=time1-FLOAT(MGASH)
                fact1=tcurves(itcurve)%dfact_curve(mgash)
                fact2=tcurves(itcurve)%dfact_curve(ngash)
                dfact=fact1*(1.0-time1)+time1*fact2
                dfact=dfact*ample
            end if                                 !hxl

        elseif(type_curve=='DEXPONENTIAL')  then   !d(f)=a*b*exp(b*t)
            time1=tcurves(itcurve)%ttime_curve(1)    !! b
            fact1=tcurves(itcurve)%dfact_curve(1)    !! a
            dfact=fact1*time1*exp(time1*time)
        elseif(type_curve=='DABT')  then  !d(f)=a*b/(b+t)**2
            time1=tcurves(itcurve)%ttime_curve(1)    !! a
            fact1=tcurves(itcurve)%dfact_curve(1)    !! b
            dfact=fact1*time1/(fact1+time)**2
        else if(type_curve=='ARCLENGTH') then
            goto 1
        else if(type_curve=='DISCONTROL') then
            goto 1
        else if(type_curve=='EXTRAPOLATION') then
            goto 1
        else
            print *,'itcurve=',itcurve, 'no type_curve=',type_curve
            stop 'no type_curve'
        end if
        tcurves(itcurve)%dfact=dfact
1       continue
    end do

    end subroutine dfact_time_curve

    subroutine read_permanent_strain !20231010
    real(irk), allocatable::sigma0(:),dmatx(:,:),strain0(:)
    character(30)text,SPtype*10
    integer(ink) ingroup,ilgroup,ielem,igaus,ie,ig,index,order_int,ngaus,nstre,nelgroup
    real(irk) mu,Emoduls
    SPtype='PE'

    read(stnunit,*)text
    if(iblks==1)element(:)%icper=0
    read(stnunit,*)ingroup

    do ilgroup=1,ingroup
        read(stnunit,*)igroup
        index=group(igroup)%index
        order_int=elkn(index)%el_field(1)%order_intrules(1)
        ngaus=elkn(index)%ggaus(order_int)%ngaus
        nstre=group(igroup)%nstre
        allocate(sigma0(nstre),dmatx(nstre,nstre),strain0(nstre))
        strain0=0.0;sigma0=0.0;dmatx=0.0
        do ie=1,group(igroup)%nelgroup
            ielem=group(igroup)%list(ie)
            if(iblks==1)then
                allocate(element(ielem)%strainx0(nstre,ngaus))
                element(ielem)%strainx0=0.0
            endif

            element(ielem)%icper=1
            do igaus=1,ngaus
                read(stnunit,*)i0,ig,strain0
                strain0=-1*strain0
                element(ielem)%strainx0(:,igaus)=strain0
            enddo	 !end do igaus
        enddo   !end do ie
        deallocate(sigma0,dmatx,strain0)
    enddo  !end do ilgroup

    end subroutine read_permanent_strain

    subroutine read_initial

    character(30)text,name,material,field1,model
    integer(ink) iinit,ingroup,igroup,ielem,igaus,i0,ie,ipoin,i1,i2,jdimn
    integer(ink) index,order_int,ngaus,nstre,nelgroup,ilgroup,matno,ngvar
    integer(ink) igaps,ipairs,npairs,unitread,kinit_g !20211214
    real(irk)    coef1,coef2,aera
    real(irk),   allocatable::sigma0(:),stres_poin(:,:),rr0(:,:),rr(:,:),  &
        tt(:,:),sgtot(:),sgloc(:),tti(:,:),sigma1(:),trot(:,:)  !20200330

    integer(ink) indofix,idofix,jdofix,totvi,ordert,vdimn,jndex,order_jnt
    real   (irk) rdofix
    real(irk),pointer::rotation(:,:)

    integer(ink) icdofn,inpoin,idofn,jpoin,idimn
    integer(ink),allocatable::ilmdofn(:),ilcdofn(:),order(:)
    real   (irk),allocatable::vinit(:)

    unitread=initunit  !20210207
    if(winit==1)unitread=initwunit !20210207
    rewind(unitread)  !20230708


    do iinit=1,ninit

        read(unitread,*)text
        write(7,*)'iinit=',iinit,'text=',text

        select case(text)

        case('STRESS')      !! for initial stresses

            read(unitread,*)ingroup,kinit_g !20211214
            do ilgroup=1,ingroup
                !if (kinit==1) then
                !   read(unitread,*)igroup,vdimn,coef1,coef2
                !else if(kinit==2) then
                read(unitread,*)igroup,vdimn,coef1,coef2
                !end if
                index=group(igroup)%index
                group(igroup)%kinit_g=kinit_g  !20211214

                write(7,*)'read_initial igroup=',igroup,'kinit_g=',kinit_g
                order_int=elkn(index)%el_field(1)%order_intrules(1)

                if(index==20.or.index==21)then
                    ngaus=1
                    ngvar=3
                    if(ndimn==3)then
                        ngvar=6
                    endif
                    nstre=ngvar
                else
                    ngaus=elkn(index)%ggaus(order_int)%ngaus

                    field1= group(igroup)%fieldid(1:1)
                    if(field1=='U')then
                        matno = group(igroup)%matno
                        material=props(matno)%mechanical%solid%material

                        if (material=='GOODMAN') then   !! 20210207
                            jndex=1
                            if (ndimn==3.and.index==9)jndex=5
                            if (ndimn==3.and.index==23)jndex=3
                            order_jnt=elkn(jndex)%el_field(1)%order_intrules(1)
                            ngaus=elkn(jndex)%ggaus(order_jnt)%ngaus
                            model=props(matno)%mechanical%solid%Goodman%model
                        endif               !! 20210207
                    endif

                    nstre=group(igroup)%nstre
                    ngvar=group(igroup)%ngvar
                endif
                matno = group(igroup)%matno
                name=props(matno)%name
                allocate(sigma0(ngvar))
                if(index==20.or.index==21)allocate(sigma1(ngvar),trot(ngvar,ngvar))
                nelgroup=group(igroup)%nelgroup
                do ielem=1,nelgroup
                    ie = group(igroup)%list(ielem)
                    do igaus=1,ngaus
                        sigma0=0.
                        !if (name=='CONTACT')then  !20210207
                        !   sigma0(1:ndimn)=.02
                        !else
                        !read(unitread,*)i0,sigma0  !20220712
                        !read(unitread,*)i0,i1,sigma0(1:nstre)  !20231215YL
                        read(unitread,*)i0,i1,sigma0(1:ngvar)  !20231215YL
                        !endif  !20210207

                        if (vdimn/=0) then !!!!!!!!!!

                            if (ndimn==2) then
                                if(vdimn==2)then
                                    sigma0(1)=coef1*sigma0(2)
                                    sigma0(4)=coef2*sigma0(2)
                                elseif(vdimn==1)then
                                    sigma0(2)=coef1*sigma0(1)
                                    sigma0(4)=coef2*sigma0(1)
                                endif
                                sigma0(3)=0.0
                            else if(ndimn==3) then
                                if(vdimn==3)then
                                    sigma0(1)=coef1*sigma0(3)
                                    sigma0(2)=coef2*sigma0(3)
                                elseif(vdimn==2)then
                                    sigma0(1)=coef1*sigma0(2)
                                    sigma0(3)=coef2*sigma0(2)
                                elseif(vdimn==1)then
                                    sigma0(2)=coef1*sigma0(1)
                                    sigma0(3)=coef2*sigma0(1)
                                endif

                                sigma0(4:6)=0.0
                            endif

                        endif
                        !write(7,*) 'ig=',igroup,'ie=',ie,'size1=',size(element(ie)%stres0(:,igaus)),'nstre=',nstre
                        !write(7,*)'sigma0=',sigma0
                        if(kinit_g==2.and.index/=20.and.index/=21)element(ie)%stres0(:,igaus)=sigma0(1:nstre)

                        if(index==20.or.index==21)then
                            rotation=>element(ie)%rotation
                            trot=0.
                            if (ndimn==2)then
                                trot(1:ndimn,1:ndimn)=transpose(rotation)
                                trot(3,3)=1.
                            else if(ndimn==3) then
                                trot(1:3,1:3)=transpose(rotation); trot(4:6,4:6)=transpose(rotation)
                            end if
                            sigma1=trot.x.sigma0
                            element(ie)%field(1)%gpvar(1:ngvar,igaus)=-sigma1
                            element(ie)%field(1)%gpvar(ngvar+1:ngvar*2,igaus)=sigma1
                            element(ie)%field(1)%gpvar0=element(ie)%field(1)%gpvar
                            if(kinit_g==2) &
                                element(ie)%stres0=element(ie)%field(1)%gpvar  !20201203
                            nullify(rotation)
                        else

                            element(ie)%field(1)%gpvar(1:ngvar,igaus)=sigma0
                            element(ie)%field(1)%gpvar0(1:ngvar,igaus)=sigma0
                        endif
                    end do
                end do
                if (model(1:3)=='FCM')then
                    do ielem=1,nelgroup
                        ie = group(igroup)%list(ielem)
                        read(unitread,*)i0,element(ie)%field(1)%strain
                        element(ie)%field(1)%strain0=element(ie)%field(1)%strain
                    enddo
                endif
                deallocate(sigma0)
                if(index==20.or.index==21)deallocate(sigma1,trot)
            end do          !!   for do ingroup

        case('INTERNAL_FORCE_BEAM') !20210207

            read(unitread,*)ingroup,kinit_g !20211214
            print *, 'ingroup=',ingroup
            DO ilgroup =1,ingroup
                read(unitread,*)igroup
                index=group(igroup)%index
                print*,'igroup=',igroup,'index=',index
                if(index/=20.and.index/=21) then
                    print *,'stop in read_initial_INTERNAL_FORCE_BEAM'
                    stop
                endif
                ngaus=1
                nstre=6
                if(ndimn==3)nstre=12
                print *,'nstre=',nstre
                ! loop for 1:nelgroup
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    read(unitread,*)i0,element(ielem)%field(1)%gpvar(1:nstre,1)
                    element(ielem)%field(1)%gpvar0=element(ielem)%field(1)%gpvar
                    if(kinit_g==2) &
                        element(ielem)%stres0=element(ielem)%field(1)%gpvar  !20201203
                end do
            end do
        case('STRESS_BOND_SLIP') !20210207

            read(unitread,*)ingroup,kinit_g !20211214
            DO ilgroup =1,ingroup
                read(unitread,*)igroup
                index=group(igroup)%index
                if(index/=25) then
                    print *,'stop in read_initial-STRESS_BOND_SLIP'
                    stop
                endif
                ngvar=2*ndimn
                ngaus=1
                nstre=ngvar
                ngvar=ngvar+5 !ngvar+1--for gaptao, ngvar+2--for gapnorm ,ngvar+3--for Ks, ngvar+4--for steel strain
                !ngvar+5--for state 0-close 1-open , integer it first! for lhg ngvar+6 ic_yty

                ! loop for 1:nelgroup
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    read(unitread,*)i0,element(ielem)%field(1)%gpvar(1:ngvar,1)
                    element(ielem)%field(1)%gpvar0=element(ielem)%field(1)%gpvar
                    if(kinit_g==2) &
                        element(ielem)%stres0(1:nstre,1)=element(ielem)%field(1)%gpvar(1:nstre,1)  !20210207
                end do

            end do

        case('CONTACT_STATE') !zhao 05/07/19 !contact
            read(unitread,*)ingroup
            DO ilgroup =1,ingroup
                read(unitread,*)igroup
                field1= group(igroup)%fieldid
                index = group(igroup)%index
                if (appear_process(igroup,iblks)>0.and.field1=='U')  then
                    matno = group(igroup)%matno
                    name  = props(matno)%name
                    if (name/='CONTACT')then
                        print *,'stop in read_initial-CONTACT'
                        stop
                    endif

                    material=props(matno)%mechanical%solid%material
                    if(material=='GOODMAN')then
                        model=props(matno)%mechanical%solid%Goodman%model
                    endif

                    !! contact
                    DO ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        read(unitread,*)i0,element(ielem)%field(1)%gapn
                        read(unitread,*)i0,element(ielem)%field(1)%gapg
                        read(unitread,*)i0,element(ielem)%field(1)%state

                        if (model(1:3)=='FCM')then
                            read(unitread,*)i0,element(ielem)%field(1)%strain
                            element(ielem)%field(1)%strain0=element(ielem)%field(1)%strain
                        endif

                        element(ielem)%field(1)%gapn0=element(ielem)%field(1)%gapn
                        element(ielem)%field(1)%gapg0=element(ielem)%field(1)%gapg
                        element(ielem)%field(1)%state0=element(ielem)%field(1)%state
                        element(ielem)%field(1)%state1=element(ielem)%field(1)%state


                        !write(chkunit,'(i10,30e14.5)')ielem,element(ielem)%field(1)%gapn
                        !write(chkunit,'(i10,30e14.5)')ielem,element(ielem)%field(1)%gapg
                        !write(chkunit,'(i10,5x,30a10)')ielem,element(ielem)%field(1)%state
                    end do
                endif
            enddo

        case('CONTACTCTT') !ctt2005 zhao 05/09/07

            if(block_stab==2)then   !20200330
                allocate(stres_poin(3*(ndimn-1),npoin))
                read(unitread,*)text
                read(unitread,*)coef1
                do ipoin=1,npoin
                    read(unitread,*)i0,stres_poin(:,ipoin)
                end do
                stres_poin=stres_poin*coef1

                allocate(rr0(ndimn,ndimn),rr(ndimn+1,ndimn+1),tt(3*(ndimn-1),3*(ndimn-1)),  &
                    sgtot(3*(ndimn-1)),sgloc(3*(ndimn-1)),tti(3*(ndimn-1),3*(ndimn-1)))

                do igaps=1,ngaps
                    npairs=gaps(igaps)%npairs
                    do ipairs=1,npairs
                        i1=gaps(igaps)%pairnode(1,ipairs)
                        i2=gaps(igaps)%pairnode(2,ipairs)
                        sgtot=.5*(stres_poin(:,i1)+stres_poin(:,i2))
                        aera=gaps(igaps)%aera(ipairs)
                        rr0=gaps(igaps)%rot(:,:,ipairs)

                        rr(1:ndimn,1:ndimn)=rr0
                        rr(1:ndimn,ndimn+1)=rr(1:ndimn,1)
                        rr(ndimn+1,1:ndimn)=rr(1,1:ndimn)
                        rr(ndimn+1,ndimn+1)=rr(1,1)

                        tt=0.0
                        tt(1:ndimn,1:ndimn)=rr**2

                        tti=0.
                        tti(1:ndimn,1:ndimn)=(transpose(rr))**2
                        if (ndimn==2) then
                            tt(1,3)=2*rr(1,1)*rr(1,2)
                            tt(2,3)=2*rr(2,1)*rr(2,2)
                            tt(3,1)=rr(1,1)*rr(2,1)
                            tt(3,2)=rr(1,2)*rr(2,2)
                            tt(3,3)=rr(1,1)*rr(2,2)+ rr(2,1)*rr(1,2)


                            tti(3,1)=rr(1,1)*rr(1,2)
                            tti(3,2)=rr(2,1)*rr(2,2)
                            tti(1,3)=2*rr(1,1)*rr(2,1)
                            tti(2,3)=2*rr(1,2)*rr(2,2)
                            tti(3,3)=rr(1,1)*rr(2,2)+ rr(2,1)*rr(1,2)
                        else if(ndimn==3) then

                            do idimn=1,ndimn
                                do jdimn=1,ndimn
                                    tt(idimn,3+jdimn)=2*rr(idimn,jdimn)*rr(idimn,jdimn+1)
                                    tt(3+idimn,jdimn)=rr(idimn,jdimn)*rr(idimn+1,jdimn)
                                    tt(3+idimn,3+jdimn)=rr(idimn,jdimn)*rr(idimn+1,jdimn+1)+  &
                                        rr(idimn+1,jdimn)*rr(idimn,jdimn+1)

                                    tti(3+jdimn,idimn)=rr(idimn,jdimn)*rr(idimn,jdimn+1)
                                    tti(jdimn,3+idimn)=2*rr(idimn,jdimn)*rr(idimn+1,jdimn)
                                    tti(3+jdimn,3+idimn)=rr(idimn,jdimn)*rr(idimn+1,jdimn+1)+  &
                                        rr(idimn+1,jdimn)*rr(idimn,jdimn+1)
                                end do
                            end do
                        endif

                        sgloc=tt.x.sgtot
                        gaps(igaps)%ctforce0(ndimn,ipairs)=sgloc(ndimn)*aera
                        if(ndimn==2)then
                            gaps(igaps)%ctforce0(1,ipairs)=sgloc(ndimn+1)*aera
                        elseif(ndimn==3)then
                            gaps(igaps)%ctforce0(1,ipairs)=sgloc(6)*aera
                            gaps(igaps)%ctforce0(2,ipairs)=sgloc(5)*aera
                        endif

                        gaps(igaps)%ctforce(:,ipairs)=gaps(igaps)%ctforce0(:,ipairs)
                        gaps(igaps)%state(ipairs)=1
                        gaps(igaps)%state0(ipairs)=1
                        !write(7,*)'ipairs=',ipairs,'ctforce=',gaps(igaps)%ctforce(:,ipairs)


                    end do
                end do
                deallocate(rr0,rr,tt,sgtot,sgloc,tti)


            else !20200330
                do igaps=1,ngaps
                    npairs=gaps(igaps)%npairs
                    do ipairs=1,npairs
                        if(kstab/=0.)then
                            read(unitread,*)i0,i0,gaps(igaps)%ctforce0(:,ipairs),gaps(igaps)%gap0(ndimn,ipairs), &
                                gaps(igaps)%state0(ipairs) !,gaps(igaps)%ft(ipairs),(gaps(igaps)%kxyz(idimn,idimn,ipairs),idimn=1,3*(ndimn-1))
                            gaps(igaps)%state0(ipairs)=1
                        else
                            read(unitread,*)i0,i0,gaps(igaps)%ctforce0(:,ipairs),gaps(igaps)%gap0(ndimn,ipairs), &
                                gaps(igaps)%state0(ipairs),gaps(igaps)%ft(ipairs),(gaps(igaps)%kxyz(idimn,idimn,ipairs),idimn=1,ndimn)
                        endif
                        gaps(igaps)%gap(ndimn,ipairs)            =gaps(igaps)%gap0(ndimn,ipairs)
                        gaps(igaps)%state(ipairs)          =gaps(igaps)%state0(ipairs)
                        gaps(igaps)%ctforce(1:ndimn,ipairs)=gaps(igaps)%ctforce0(1:ndimn,ipairs)
                        !write(7,*)'ipairs=',ipairs,'ctforce=', gaps(igaps)%ctforce(1:ndimn,ipairs),'state=',gaps(igaps)%state(ipairs)
                        if(kinit==2)gaps(igaps)%ctforce_stres0(1:ndimn,ipairs)=gaps(igaps)%ctforce0(1:ndimn,ipairs) !2019/03/19
                        !if(gaps(igaps)%state(ipairs)==2.and.xlwsol==1)then
                        ! gaps(igaps)%kxyz(1,1,ipairs)=gaps(igaps)%kgroup1(1,1)  !柔度系数取大值模拟自由滑动
                        !if(ndimn==3)gaps(igaps)%kxyz(2,2,ipairs)=gaps(igaps)%kgroup1(2,2)  !柔度系数取大值模拟自由滑动
                        !end if
                    enddo
                enddo
            endif !20200330
        case('REACTION')      !! for initial reactions

            read(unitread,*)indofix
            do jdofix=1,indofix
                read(unitread,*)idofix,rdofix
                prescrib(idofix)%rdofix=rdofix
            end do

            case default    !! for initial(Ux,Uy,Uz,Thxy,....)
            !! it could be velocity or acceleration
            allocate(ilmdofn(mdofn),ilcdofn(mdofn),order(1:mdofn))
            read(unitread,*)i0,inpoin
            read(unitread,*)ilmdofn(1:mdofn)
            read(unitread,*)order(1:mdofn)
            icdofn=0
            do idofn=1,mdofn
                if (ilmdofn(idofn)/=0) then
                    icdofn=icdofn+1
                    ilcdofn(icdofn)=idofn
                endif
            end do
            allocate(vinit(icdofn))

            !write(7,*)'ipoin,idofn,totvi,result_zero(totvi)='
            do jpoin=1,inpoin
                read(unitread,*)ipoin,vinit

                do idofn=1,icdofn
                    totvi=nodfn(lmdofn(ilcdofn(idofn)),ipoin)

                    if(ilcdofn(idofn)==8) then  !20230331
                        if(vinit(idofn)<=0)vinit(idofn)=0.
                    endif !20230331
                    ordert=order(ilcdofn(idofn))
                    if (ordert==0.and.totvi/=0)then
                        result_zero(totvi)=vinit(idofn)
                        !write(7,*)ipoin,idofn,totvi,result_zero(totvi)

                    elseif(ordert==1)then
                        result_first(totvi)=vinit(idofn)
                    elseif(ordert==2)then
                        result_second(totvi)=vinit(idofn)
                    endif

                    !            if(ilcdofn(idofn)/=8.and.ordert==0.and.totvi/=0) result_zero(totvi)=0. !special zhao

                    !! initial pore pressure
                    !if(type_problem=='F'.and.(ilcdofn(idofn)==8.and.ordert==0)) prstat(ipoin)=vinit(idofn) !20221013
                    if(ilcdofn(idofn)==8.and.ordert==0) prstat(ipoin)=vinit(idofn) !20221013
                end do
            end do
1           continue
            deallocate(vinit,order)

        end select

    end do     !!!iinit

    end subroutine read_initial

    subroutine placement_temperature(result0)
    character(10) fieldid
    integer(ink) igroup,ic,ifield,matno,place_curve,ielgroup,ipoin,jpoin,nline_g_w,npairs_wc,  &
        ielem, nrfields,index,nnode,itotv,jtotv,iintf,nintf,i0,i1,iwc,twater_curve
    integer(ink), pointer::ldofs(:),pairnode_wc(:)
    real(irk),allocatable::result1(:)
    real(irk) ::result0(ntotv),dfact

    real   (irk) temperature,tp,tpave,tpele


    allocate(result1(ntotv))
    result0=result_zero;
    result1=result_zero;
    print*,'maxval(result_zero)=',maxval(result_zero)
    DO igroup =1,ngroup
        if(appear(igroup)<=0)cycle
        if(appear_process(igroup,iblks-1)==0.and.   &
            appear_process(igroup,iblks)==1) then

            nrfields=group(igroup)%nrfields
            fieldid=group(igroup)%fieldid

            ic=0

            do ifield=1,nrfields
                if(fieldid(ifield:ifield)=='T')then
                    ic=1
                    exit
                end if
            end do

            if(ic==1) then
                matno = group(igroup)%matno
                place_curve=props(matno)%heat%place_curve
                if(place_curve/=0) then
                    temperature=tcurves(place_curve)%dfact
                    group(igroup)%water_pipe%tp=temperature
                    write(chkunit,*)'temp=',temperature
                    DO ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        ldofs => element(ielem)%field(ifield)%ldofs_f
                        result_zero(ldofs)=temperature

                        !do idofn=1,size(ldofs)     !修改浇筑层界面温度
                        !	if(abs(result_zero(ldofs(idofn)))<.0001) then
                        !		result_zero(ldofs(idofn))=temperature
                        !	else
                        !		result_zero(ldofs(idofn))=(result0(ldofs(idofn))+temperature)*.5
                        !	endif
                        !            enddo

                        do idofn=1,size(ldofs)     !cj042@126.com  20191120修改浇筑层界面温度
                            if(abs(result0(ldofs(idofn)))<.0001) then
                                result_zero(ldofs(idofn))=temperature
                                result1(ldofs(idofn))=temperature  !20191127修改浇筑层除接触面意外的节点的基准温度
                            else
                                result_zero(ldofs(idofn))=(result0(ldofs(idofn))+temperature)*.5  !基准值result0不变【以下层为准】
                            endif
                        enddo

                        nullify(ldofs)
                    end do

                    index=group(igroup)%index
                    nnode=elkn(index)%el_field(ifield)%nnode_f
                    tpave=0.
                    DO ielgroup = 1,group(igroup)%nelgroup      !ielgroup
                        ielem = group(igroup)%list(ielgroup)
                        ldofs =>element(ielem)%field(1)%ldofs_f
                        tpele=sum(result_zero(ldofs))/nnode
                        tpave=tpave+tpele
                        nullify(ldofs)
                    end do
                    tp=tpave/group(igroup)%nelgroup
                    group(igroup)%water_pipe%tp=tp
                    group(igroup)%water_pipe%time_eq=0  ! cj042 20191104 equivalent age
                    group(igroup)%water_pipe%time_real=0  ! cj042 20191104 equivalent age

                endif ! endif for place_curve/=0
            endif ! endif for any fileld=='T'(ic=1)
        endif ! endif for the group just appears in the block
    end do ! igroup
    result0=result1
    deallocate(result1)

    !! 对从节点数值进行修改，以保证水管单元初始温度与依赖的混凝土节点温度相同(20210417)
    do itotv=1,ntotv
        nintf=trans(itotv)%nintf
        if (nintf==0) cycle
        result_zero(itotv)=0.
        do iintf=1,nintf
            jtotv=trans(itotv)%listf(iintf)
            result_zero(itotv)=result_zero(itotv)+result_zero(jtotv)*trans(itotv)%rintf(iintf)
        end do
    end do

    do i0=1,nwcpipe
        iwc=wc_pipe(i0)%iwc
        twater_curve=wc_pipe(i0)%twater_curve
        dfact   =tcurves(twater_curve)%dfact
        nline_g_w=wc_pipe(i0)%nline_g_w
        do i1=1,nline_g_w
            npairs_wc=wc_pipe(i0)%line_g_w(i1)%npairs_wc
            pairnode_wc=>wc_pipe(i0)%line_g_w(i1)%pairnode_wc
            if(iwc==1)then
                ipoin=pairnode_wc(1)
                jtotv=nodfn(lmdofn(10),ipoin)
                result_zero(jtotv)=dfact
            elseif(iwc==-1)then
                ipoin=pairnode_wc(npairs_wc)
                jtotv=nodfn(lmdofn(10),ipoin)
                result_zero(jtotv)=dfact
            endif
            nullify(pairnode_wc)
        end do
    end do
    !!end 对从节点数值修改，以保证水管单元初始温度与依赖的混凝土节点温度相同(20210417)



    end   subroutine placement_temperature

    subroutine placement_temperature0

    character(10) fieldid
    integer(ink) igroup,ic,ifield,matno,place_curve,ielgroup,ielem, nrfields
    integer(ink), pointer::ldofs(:)
    real   (irk) temperature
    DO igroup =1,ngroup
        if(appear(igroup)<=0)cycle
        if (appear_process(igroup,iblks-1)==0.and.   &
            appear_process(igroup,iblks)==1) then
            nrfields=group(igroup)%nrfields
            fieldid=group(igroup)%fieldid
            ic=0
            do ifield=1,nrfields
                if (fieldid(ifield:ifield)=='T')then
                    ic=1
                    exit
                end if
            end do
            if (ic==1) then
                matno = group(igroup)%matno
                place_curve=props(matno)%heat%place_curve
                if (place_curve/=0) then
                    temperature=tcurves(place_curve)%dfact
                    DO ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        ldofs => element(ielem)%field(ifield)%ldofs_f
                        result_zero(ldofs)=temperature
                        nullify(ldofs)
                    end do
                endif ! endif for place_curve/=0
            endif ! endif for any fileld=='T'(ic=1)
        endif ! endif for the group just appears in the block
    end do ! igroup

    end   subroutine placement_temperature0

    SUBROUTINE modf_time_order
    character(10)fieldid,name
    integer(ink) index,matno

    do igroup=1,ngroup
        fieldid=group(igroup)%fieldid
        index=group(igroup)%index
        matno = group(igroup)%matno

        !print *,'igroup=',igroup   !20231215YL

        if(fieldid=='U') then
            elkn(index)%el_field(1)%order_time(1)=0
            elkn(index)%el_field(1)%order_time(2)=2
            if(type_problem=='S')elkn(index)%el_field(1)%order_time(2)=1
        elseif(fieldid=='W') then
            elkn(index)%el_field(1)%order_time(1)=0
            elkn(index)%el_field(1)%order_time(2)=2
            name=props(matno)%name
            if(name(1:6)=='NSSoil')elkn(index)%el_field(1)%order_time(2)=1

        endif

        if(fieldid=='T'.and.type_problem=='S') then
            elkn(index)%el_field(1)%order_time(1)=0
            elkn(index)%el_field(1)%order_time(2)=1
        endif

        if(fieldid=='UW') then
            elkn(index)%el_field(2)%order_time(1)=0
            elkn(index)%el_field(2)%order_time(2)=1


            elkn(index)%couple(1)%order_couple(1)=0
            elkn(index)%couple(1)%order_couple(2)=0

            if(type_problem/='Q') then
                elkn(index)%couple(1)%order_couple(1)=1
                elkn(index)%couple(1)%order_couple(2)=0
            endif

            if(type_problem=='S') then
                elkn(index)%el_field(1)%order_time(1)=0
                elkn(index)%el_field(1)%order_time(2)=1

            elseif(type_problem=='F') then
                elkn(index)%el_field(1)%order_time(1)=0
                elkn(index)%el_field(1)%order_time(2)=2
            endif

        endif

    end do

    END  SUBROUTINE modf_time_order
