    SUBROUTINE unit_force_trans(igapb,kdimn,ij,unitg,igaps,ipairs,rvector)
    real(irk) coef1,rvector(:),unitg(:)
    integer(ink) igapb,igaps,ipairs,ij1,ij2,j0,ij,jpoin,jdimn,jtotv,nintf,iieq,iintf,kdimn,nnodei,nnodej

    nnodei=size(gaps(igaps)%pairnode(:,ipairs))  !2017/02/14
    nnodej=nnodei/2    !2017/02/14
    if(ij==1)then
        ij1=1
        ij2=1
        if(contactpe==2)then
            ij1=1
            !ij2=2*(ndimn-1)
            ij2=nnodej    !2017/02/14
        endif
    elseif(ij==2)then
        ij1=2
        ij2=2
        if(contactpe==2)then
            !ij1=2*(ndimn-1)+1
            !ij2=4*(ndimn-1)
            ij1=nnodej+1   !2017/02/14
            ij2=nnodei     !2017/02/14
        endif
    endif

    coef1=1.
    !if(contactpe==2)coef1=.5/(ndimn-1)
    if(contactpe==2)coef1=1./nnodej  !2017/02/14
    do j0=ij1,ij2
        jpoin=gaps(igaps)%pairnode(j0,ipairs)


        do jdimn=1,kdimn
            jtotv=nodfn(jdimn,jpoin)
            nintf=trans(jtotv)%nintf
            if(nintf/=0) then
                iieq=totveq(jtotv)
                if(iieq/=0)rvector(iieq)=0.
                do iintf=1,nintf
                    iieq=totveq(trans(jtotv)%listf(iintf))
                    if(iieq/=0) &
                        rvector(iieq)=rvector(iieq)+coef1*unitg(jdimn)*trans(jtotv)%rintf(iintf)
                end do
            else
                if(totveq(jtotv)/=0) &
                    rvector(totveq(jtotv))=rvector(totveq(jtotv))+coef1*unitg(jdimn)

            endif
        end do !jdimn
    end do
    end SUBROUTINE unit_force_trans

    SUBROUTINE unit_dis_force_trans(igapb,jpoin,kdimn,unitg,rvector) !20210820

    real(irk) rvector(:),unitg(:),dispoint
    real(irk),allocatable::dis_unit(:),load_unit(:),eload(:),value(:)
    integer(ink) kdimn,jpoin,igapb,igroup,jgroup,idofn,itotv,nintf,iieq,iintf,nevab, &
        ielem,ielgroup,jtotv
    integer(ink),pointer::ldofs(:)
    real(irk),   pointer::fstif(:,:)

    allocate(dis_unit(ntotv),load_unit(ntotv))
    dis_unit=0.;load_unit=0.

    do idofn=1,kdimn
        itotv=nodfn(idofn,jpoin)
        dis_unit(itotv)=unitg(idofn)
    end do

    do itotv=1,ntotv
        nintf=trans(itotv)%nintf
        if (nintf==0) cycle
        dispoint=0.
        do iintf=1,nintf
            jtotv=trans(itotv)%listf(iintf)
            dispoint=dispoint+dis_unit(jtotv)*trans(itotv)%rintf(iintf)
        end do
        dis_unit(itotv)=dispoint
    end do


    do igroup=1,gapb(igapb)%ngroupb
        jgroup=gapb(igapb)%listgroupb(igroup)
        if (appear(jgroup)<=0) cycle
        nevab=size(element(group(jgroup)%list(1))%field(1)%ldofs_f)
        allocate(value(nevab),eload(nevab))
        nevab=0.;eload=0.
        DO ielgroup = 1,group(jgroup)%nelgroup
            ielem = group(jgroup)%list(ielgroup)

            if (associated(element(ielem)%field(1)%khandmc(1)%fstif)) then
                fstif=>element(ielem)%field(1)%khandmc(1)%fstif
                ldofs=>element(ielem)%field(1)%ldofs_f
                value=dis_unit(ldofs)
                eload=fstif.x.value
                load_unit(ldofs)=load_unit(ldofs)-eload
                nullify(fstif,ldofs)
            endif
        end do       !!ielgroup
        deallocate(value,eload)
    end do     !!  for igroup

    do itotv=1,ntotv
        if (totveq(itotv)==0)cycle
        rvector(totveq(itotv))=rvector(totveq(itotv))+load_unit(itotv)
    end do

    !!int2000
    do itotv=1,ntotv
        nintf=trans(itotv)%nintf
        if (nintf==0) cycle
        do iintf=1,nintf
            iieq=totveq(trans(itotv)%listf(iintf))
            if(iieq/=0)rvector(iieq)=rvector(iieq)+  &
                load_unit(itotv)*trans(itotv)%rintf(iintf)
        end do
    end do
    !!int2000
    deallocate(dis_unit,load_unit)
    end SUBROUTINE unit_dis_force_trans  !20210820

    subroutine result_node_to_center(kdimn,ij,igaps,ipairs,dx,dis)
    integer(ink) igaps,ipairs,ij1,ij2,j0,i12,ij,kdimn,nnodei,nnodej
    real(irk) coef1
    real(irk) dx(:),dis(:)

    nnodei=size(gaps(igaps)%pairnode(:,ipairs))  !2017/02/14
    nnodej=nnodei/2    !2017/02/14
    if(ij==1)then
        ij1=1
        ij2=1
        if(contactpe==2)then
            ij1=1
            !ij2=2*(ndimn-1)
            ij2=nnodej    !2017/02/14
        endif
    elseif(ij==2)then
        ij1=2
        ij2=2
        if(contactpe==2)then
            !ij1=2*(ndimn-1)+1
            !ij2=4*(ndimn-1)
            ij1=nnodej+1   !2017/02/14
            ij2=nnodei     !2017/02/14
        endif
    endif

    coef1=1.
    !if(contactpe==2)coef1=.5/(ndimn-1)
    if(contactpe==2)coef1=1./nnodej  !2017/02/14


    dis=0.
    do j0=ij1,ij2
        i12=gaps(igaps)%pairnode(j0,ipairs)
        dis=dis+coef1*dx(nodfn(1:kdimn,i12))
    end do  !j0
    end subroutine  result_node_to_center

    subroutine contact_state(ic)
    character(20) field1,name,material
    integer (ink) igroup,matno,index,nnode,nnode_half,  &
        order_int,ngaus,ielgroup,ielem,inode, &
        igaus,nevab,idofn,igap0,jndex,ic,iiii,jjjj,icftg,icft
    integer (ink) gap_kind,ngapx,ix  !20231006
    real    (irk) eps,gap0,smean,ftcontact,ft0,radius,alfax,dxi  !20231006
    integer (ink),pointer::ldofs(:)
    real    (irk),allocatable::rot(:),eldis(:),nordis(:),shape(:,:),gapnod(:),gapgaus(:)
    real    (irk),allocatable::centerx(:),gapx(:),gapalfax(:)  !20231006
    real    (irk):: current_gap, gap_change, natural_gap
    character(20):: previous_state, current_state
    real    (irk):: old_gap


    eps = 1.e-5  ! 增大eps提高数值稳定性


    write(7,*)'contact_state********************'
    DO igroup =1,ngroup
        field1= group(igroup)%fieldid(1:1)
        if (appear(igroup)>0.and.field1=='U')then
            matno = group(igroup)%matno
            name=props(matno)%name
            material=props(matno)%mechanical%solid%material
            igap0=0
            if (name=='CONTACT')then
                gap0     =props(matno)%mechanical%solid%gap0
                igap0    =props(matno)%mechanical%solid%igap0

                if(igap0==2) then    !20231006
                    ngapx=props(matno)%mechanical%solid%gap_define%ngapx
                    gap_kind=props(matno)%mechanical%solid%gap_define%gap_kind
                    if(gap_kind==1) then
                        allocate(centerx(ndimn),gapx(ngapx),gapalfax(ngapx))
                        radius=props(matno)%mechanical%solid%gap_define%radius
                        centerx=props(matno)%mechanical%solid%gap_define%centerx
                        gapalfax=props(matno)%mechanical%solid%gap_define%gapalfax
                        gapx=props(matno)%mechanical%solid%gap_define%gapx
                    endif
                endif !20231006

                ftcontact=props(matno)%mechanical%solid%ftcontact !zhao 05/08/02
                icft     =props(matno)%mechanical%solid%icft
                index    =group(igroup)%index
                nnode    =elkn(index)%el_field(1)%nnode_f
                nnode_half=nnode/2
                order_int=elkn(index)%el_field(1)%order_intrules(1)
                ngaus    =elkn(index)%ggaus(order_int)%ngaus
                if(material=='GOODMAN') then   !! for goodman element

                    jndex=1
                    if (ndimn==3.and.index==9)jndex=5  !2017/02/14
                    if (ndimn==3.and.index==23)jndex=3  !2017/02/14
                    order_int=elkn(jndex)%el_field(1)%order_intrules(1)
                    ngaus=elkn(jndex)%ggaus(order_int)%ngaus
                endif !! for goodman element
                nevab=nnode*ndimn
                allocate(eldis(nevab),nordis(nnode),gapnod(nnode),gapgaus(ngaus),rot(ndimn))
                if (material=='GOODMAN') then   !! for goodman element
                    allocate(shape(nnode_half,ngaus))
                else
                    allocate(shape(nnode,ngaus))
                endif

                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    if (tension_contact(ielem)/=1)then !zhao 05/07/22  tcl,original==1
                        if (material=='GOODMAN') then
                            rot=element(ielem)%rotation(ndimn,:)
                        else
                            rot=element(ielem)%rotation(1,:)
                        endif
                        ldofs => element(ielem)%field(1)%ldofs_f
                        !eldis=delitfi(ldofs)
                        eldis = deltafi(ldofs)
                        !                  eldis = result_zero(ldofs)

                        if(appear_process(igroup,iblks)==1.and.     &
                            appear_process(igroup,iblks-1)==0.and.   &
                            iincs==1.and.istep==inc_step.and.iiter==1.and.ic==0) then

                            if (igap0==1) then
                                element(ielem)%field(1)%gapg =gap0
                                element(ielem)%field(1)%gapg0=gap0
                                element(ielem)%field(1)%gapn =gap0
                                element(ielem)%field(1)%gapn0=gap0
                            elseif(igap0==2.and.gap_kind==1)then !20231006
                                gapnod=0.
                                do inode=1,nnode
                                    idofn=element(ielem)%field(1)%lnods_f(inode)
                                    dxi=coord(1,idofn)-centerx(1)
                                    alfax=dxi/radius
                                    alfax=acosd(alfax)
                                    if(alfax<gapalfax(1).or.alfax>gapalfax(ngapx))cycle
                                    do ix=1,ngapx-1
                                        if(alfax>=gapalfax(ix).and.alfax<=gapalfax(ix+1))then
                                            gapnod(inode)=gapx(ix)+(gapx(ix+1)-gapx(ix))*(alfax-gapalfax(ix))/(gapalfax(ix+1)-gapalfax(ix))
                                        endif
                                    end do
                                end do

                                if (material=='GOODMAN') then
                                    shape = elkn(jndex)%ggaus(order_int)%shape(:,:)
                                    gapgaus=transpose(shape).x.gapnod(1:nnode_half)
                                else
                                    shape = elkn(index)%ggaus(order_int)%shape(:,:)
                                    gapgaus=transpose(shape).x.gapnod
                                endif
                                !write(7,*)'ie=',ielem,'gapgaus=',gapgaus
                                element(ielem)%field(1)%gapg0=gapgaus
                                element(ielem)%field(1)%gapg=element(ielem)%field(1)%gapg0
                                element(ielem)%field(1)%gapn=element(ielem)%field(1)%gapn0
                            elseif(igap0==99)then !zhao 05/07/19
                                do inode=1,nnode
                                    idofn=element(ielem)%field(1)%lnods_f(inode)
                                    nordis(inode)=rot.d.coord(:,idofn)
                                end do
                                if (ndimn==3) then
                                    do inode=1,nnode_half
                                        gapnod(inode)=nordis(inode+nnode_half)-nordis(inode)
                                        gapnod(inode+nnode_half)=gapnod(inode)
                                    end do
                                else  ! 2D
                                    gapnod(1)=nordis(4)-nordis(1)
                                    gapnod(4)=gapnod(1)
                                    gapnod(2)=nordis(3)-nordis(2)
                                    gapnod(3)=gapnod(2)
                                endif
                                element(ielem)%field(1)%gapn0=gapnod

                                if (material=='GOODMAN') then
                                    shape = elkn(jndex)%ggaus(order_int)%shape(:,:)
                                    gapgaus=transpose(shape).x.gapnod(1:nnode_half)
                                else
                                    shape = elkn(index)%ggaus(order_int)%shape(:,:)
                                    gapgaus=transpose(shape).x.gapnod
                                endif
                                !write(7,*)'ie=',ielem,'gapgaus=',gapgaus
                                element(ielem)%field(1)%gapg0=gapgaus
                                element(ielem)%field(1)%gapg=element(ielem)%field(1)%gapg0
                                element(ielem)%field(1)%gapn=element(ielem)%field(1)%gapn0
                                if(istep.eq.1.and.iincs==1) then
                                    element(ielem)%field(1)%natural_thickness = element(ielem)%field(1)%gapg
                                endif
                            endif  !igap0=1
                        else  !for ic==1
                            do inode=1,nnode
                                idofn=(inode-1)*ndimn
                                nordis(inode)=rot.d.eldis(idofn+1:idofn+ndimn)
                            end do
                            if (ndimn==3) then
                                do inode=1,nnode_half
                                    gapnod(inode)=nordis(inode+nnode_half)-nordis(inode)
                                    gapnod(inode+nnode_half)=gapnod(inode)
                                end do
                            else  ! 2D
                                gapnod(1)=nordis(4)-nordis(1)
                                gapnod(4)=gapnod(1)
                                gapnod(2)=nordis(3)-nordis(2)
                                gapnod(3)=gapnod(2)
                            endif

                            element(ielem)%field(1)%gapn=element(ielem)%field(1)%gapn0+gapnod

                            if (material=='GOODMAN') then
                                shape = elkn(jndex)%ggaus(order_int)%shape(:,:)
                                gapgaus=transpose(shape).x.gapnod(1:nnode_half)
                            else
                                shape = elkn(index)%ggaus(order_int)%shape(:,:)
                                gapgaus=transpose(shape).x.gapnod
                            endif
                            element(ielem)%field(1)%gapg=gapgaus+element(ielem)%field(1)%gapg0
                        endif
                        !write(7,"(A15,I10,5(A10,3E15.7))")'当前间隙为：ie=',ielem,'gapg=   ',element(ielem)%field(1)%gapg
                        !write(7,"(A15,I10,5(A10,3E15.7))")'当前间隙为：ie=',ielem,'gapg0=  ',element(ielem)%field(1)%gapg0
                        !write(7,"(A15,I10,5(A10,3E15.7))")'当前间隙为：ie=',ielem,'gapgaus=',gapgaus
                        do igaus=1,ngaus
                            icftg=element(ielem)%field(1)%icftcontact(igaus)
                            if (material=='GOODMAN') then
                                smean=element(ielem)%field(1)%gpvar(ndimn,igaus)
                            else
                                smean=element(ielem)%field(1)%ntstress(1,igaus)
                            endif
                            element(ielem)%field(1)%state1(igaus)=element(ielem)%field(1)%state(igaus) !zhao 05/07/22
                            element(ielem)%field(1)%state(igaus)='contact'
                            if (icftg==0.and.icft/=0)then
                                ft0=ftcontact
                            else
                                ft0=0.01
                            endif


                            ! 获取间隙变化量（相对于初始状态的变化）
                            current_gap = element(ielem)%field(1)%gapg(igaus)  ! 总间隙（用于显示）
                            gap_change = current_gap-element(ielem)%field(1)%gapg0(igaus)  ! 间隙变化量
                            natural_gap = element(ielem)%field(1)%natural_thickness(igaus)  ! 初始间隙


                            previous_state = element(ielem)%field(1)%state1(igaus)

                            !if (element(ielem)%field(1)%gapg(igaus)>eps.and.smean>ft0) then !ooo
                            !    element(ielem)%field(1)%icftcontact(igaus)=1
                            !    element(ielem)%field(1)%state(igaus)='open'
                            !    !                        write(7,*)'ielem=',ielem,'igaus=',igaus,'gapg=',element(ielem)%field(1)%gapg(igaus),'smean=',smean
                            !endif
                            !cycle


                            if (previous_state == 'contact') then
                                ! contact → open 判断：基于位移增量和真实应力

                                ! 修正的判断条件：
                                ! 1. 位移增量明显为正（张开方向）
                                ! 2. 应力超过拉伸极限（这里smean还是真实应力）
                                if (smean > ft0) then
                                    element(ielem)%field(1)%state(igaus) = 'open'
                                    element(ielem)%evk(:,igaus) = 0.02
                                else
                                    element(ielem)%field(1)%state(igaus) = 'contact'
                                endif

                            elseif (previous_state == 'open') then
                                ! open → contact 判断：基于间隙闭合几何条件

                                ! 正确的判断条件：
                                ! open → contact 的核心：间隙足够小，接近真正接触
                                ! 判断依据：当前间隙接近接触阈值

                                if (gap_change < eps .and. current_gap<natural_gap) then
                                    element(ielem)%field(1)%state(igaus) = 'contact'
                                else
                                    element(ielem)%field(1)%state(igaus) = 'open'
                                endif
                            endif


                            ! --- 根据最终的 'state' 字符串更新整数标志位 icftcontact ---
                            if (element(ielem)%field(1)%state(igaus) == 'open') then
                                element(ielem)%field(1)%icftcontact(igaus) = 1 ! 1 代表张开
                            else
                                element(ielem)%field(1)%icftcontact(igaus) = 0 ! 0 代表接触
                            endif

                        end do  !!igaus
                        element(ielem)%field(1)%icok=0
                        if(all(element(ielem)%field(1)%state1==element(ielem)%field(1)%state))element(ielem)%field(1)%icok=1
                        nullify(ldofs)
                    else !zhao 05/07/22
                        do igaus=1,ngaus
                            element(ielem)%field(1)%state (igaus)='open'
                            element(ielem)%field(1)%state0(igaus)='open'
                            element(ielem)%field(1)%state1(igaus)='open'
                        enddo
                    endif !zhao 05/07/22
                end do   !! ielgroup
                deallocate(rot,nordis,shape,eldis,gapnod,gapgaus)
            endif

            if (name=='CONTACT'.and.igap0==2.and.gap_kind==1)  &
                deallocate(centerx,gapx,gapalfax)   !20231006

        endif
    end do  !! igroup


    iiii=0 ; jjjj=0
    DO igroup =1,ngroup
        field1= group(igroup)%fieldid(1:1)
        if (appear(igroup)>0.and.field1=='U')then
            matno = group(igroup)%matno
            name=props(matno)%name
            material=props(matno)%mechanical%solid%material
            if (name=='CONTACT')then
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    !if(tension_contact(ielem)/=1)cycle  !20230922
                    if(tension_contact(ielem)==1)cycle  !20230922
                    jjjj=jjjj+1
                    if(element(ielem)%field(1)%icok==1)iiii=iiii+1
                enddo
            endif
        endif
    enddo
    if(iiii==jjjj)iccontact=1

    end subroutine contact_state

    subroutine crack_state

    character(20) field1,material,criteria,name
    integer(ink) iielem,iigaus,igroup,matno,index,order_int,ngaus,igaus,ielem
    real   (irk) smin,ft,sx,sy,sxy,s1,delta,strem(3),rr(3,3)
    real   (irk),allocatable::stemp(:)
    allocate(stemp(3*(ndimn-1)))
    smin=1.0e-10
    iielem=0
    iigaus=0
    DO igroup =1,ngroup
        field1= group(igroup)%fieldid(1:1)
        if (appear(igroup)>0.and.field1=='U')then
            matno = group(igroup)%matno
            name=props(matno)%name
            index = group(igroup)%index
            if (name/='CRACK'.or.index==20.or.index==21.or.index==25)cycle
            order_int=elkn(index)%el_field(1)%order_intrules(1)
            ngaus = elkn(index)%ggaus(order_int)%ngaus

            DO ielgroup = 1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                !do igaus=1,ngaus/2 !special for xld
                do igaus=1,ngaus
                    stemp=element(ielem)%field(1)%gpvar(1:3*(ndimn-1),igaus)

                    if (ndimn==3)then
                        call stresmr ( stemp, strem, rr)
                        s1=strem(1)
                    else
                        sx=stemp(1)
                        sy=stemp(2)
                        sxy=stemp(3)
                        delta=sqrt((sx-sy)**2/4+sxy**2)
                        s1=0.
                        if(delta>1.e-5)s1=(sx+sy)/2.+delta
                    endif
                    element(ielem)%field(1)%ntstress(1,igaus)=s1
                    if (smin<s1.and.s1>=ftcrack)then
                        smin=s1
                        iielem=ielem
                        iigaus=igaus
                    endif
                enddo
            enddo
        endif
    enddo
    deallocate(stemp)
    if(iielem*iigaus/=0)then
        write(7,*)'*******************************crack***************************'
        write(7,*)'istep=',istep,'  iiter=',iiter
        write(7,*)'ielem=',iielem,' igaus=',iigaus
        element(iielem)%field(1)%state(iigaus)='open'

        !if(ndimn==3)then
        !  element(iielem)%field(1)%state(iigaus+4)='open'
        !else
        !  element(iielem)%field(1)%state(iigaus+2)='open'
        !endif

        !one element failure in each step
        !element(iielem)%field(1)%state='open'

    endif

    end subroutine crack_state

    subroutine local_stress
    character(20) field1,material,criteria,name
    integer (ink) igroup,matno,index,igaus,idimn,       &
        order_int,ngaus,ielgroup,ielem,ic
    real    (irk) smean,steff
    real    (irk),allocatable::rot(:),devia(:),stemp(:),tensor(:,:)

    DO igroup =1,ngroup
        field1= group(igroup)%fieldid(1:1)
        if (appear(igroup)>0.and.field1=='U')then

            matno = group(igroup)%matno
            material=props(matno)%mechanical%solid%material
            name  = props(matno)%name
            ic=0
            if(name=='CONTACT'.and.material/='GOODMAN')ic=1
            if (material=='CLASSICALEP') then
                criteria=props(matno)%mechanical%solid%ClassicalEP%criteria
                if(criteria=='MCJOINT')ic=1
            endif
            if (ic==1) then
                allocate(rot(ndimn),devia(ndimn),stemp(3*(ndimn-1)),tensor(ndimn,ndimn))
                index = group(igroup)%index
                order_int=elkn(index)%el_field(1)%order_intrules(1)
                ngaus = elkn(index)%ggaus(order_int)%ngaus
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    rot=element(ielem)%rotation(1,:)
                    do igaus=1,ngaus
                        stemp=element(ielem)%field(1)%gpvar(1:3*(ndimn-1),igaus)
                        do idimn=1,ndimn
                            tensor(idimn,idimn)=stemp(idimn)
                        end do
                        if (ndimn==2) then
                            tensor(1,2)=stemp(3)
                            tensor(2,1)=stemp(3)
                        elseif(ndimn==3) then
                            tensor(1,2)=stemp(4)
                            tensor(1,3)=stemp(6)
                            tensor(2,3)=stemp(5)
                            tensor(2,1)=stemp(4)
                            tensor(3,1)=stemp(6)
                            tensor(3,2)=stemp(5)
                        endif
                        devia(1:ndimn)=tensor.x.rot(1:ndimn)
                        smean=rot(1:ndimn).d.devia(1:ndimn)
                        steff=sum(devia(1:ndimn)**2)-smean**2
                        if (steff.le.1.e-5) then
                            steff=1.e-5
                        else
                            steff=sqrt(steff)
                        endif
                        !                  write(chkunit,*)'smean=',smean,'steff=',steff
                        element(ielem)%field(1)%ntstress(1,igaus)=smean
                        element(ielem)%field(1)%ntstress(2,igaus)=steff
                    end do  ! igaus
                end do  ! ielgroup
                deallocate(rot,tensor,stemp,devia)
            endif    ! for ic=1
        endif  ! for group appearing
    end do  !! igroup

    end subroutine locaL_stress
