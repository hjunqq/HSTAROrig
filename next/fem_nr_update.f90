    SUBROUTINE PREDICT

    integer(ink) ordert,idofn,ipoin,itotv,idimn,jdofn
    real   (irk) accmid
    integer(ink) igroup,ielem,ielgroup,index,igapb,jdimn !simo_rifai
    character(1 )field1  !simo_rifai
    character(10)special !simo_rifa
    real(irk),allocatable::disl(:)

    !! the following is for the predicted value

    write(7,*) 'in predict'
    do idofn=1,cdofn

        jdofn=lcdofn(idofn)               !20200226
        if(jdofn==10.and.outintr/=0)cycle !20200226

        ordert=order_time_mdofn(lcdofn(idofn))


        do ipoin=1,npoin
            itotv=nodfn(idofn,ipoin)
            ! if(ipoin==168) &
            !print *,'ipoin=',ipoin,'idofn=',idofn,'itov=',itotv,'iffix=',iffix(itotv)
            if(itotv==0)cycle
            !  if(itotv/=0.and.iffix(itotv)==0)then  !!  y1
            if (itotv/=0.and.(iffix(itotv)==0.or.iffix(itotv)==4.or.iffix(itotv)==5.or.trans(itotv)%nintf/=0))then  !!  y1
                if(ordert==1)then
                    result_zero(itotv)=result_zero(itotv)+   &
                        ditime*result_first(itotv)
                    deltafi(itotv)=ditime*result_first(itotv)
                elseif(ordert==2) then
                    deltafi(itotv)=ditime*result_first(itotv)              &
                        +.5*ditime**2*result_second(itotv)
                    result_zero(itotv)=result_zero(itotv)+   &
                        ditime*result_first(itotv)+       &
                        .5*ditime**2*result_second(itotv)
                    result_first(itotv)=result_first(itotv)+ditime*result_second(itotv)
                    !if(iffix(itotv)==4)then
                    !write(7,*)'itotv=',itotv,'delta=',deltafi(itotv),'dis=',result_zero(itotv),'vec=',result_first(itotv),'acc=',result_second(itotv)
                    !endif


                endif

            else if(itotv/=0.and.(iffix(itotv)/=0.and.trans(itotv)%nintf==0))then !!! 07/04/16



                if(nbackdT==2.and.iffix(itotv)==1.and.ordert==1)then !20230216
                    result_zero(itotv)=result_zero(itotv)+result_first(itotv)*ditime
                    !write(7,*)'itotv=',itotv,'result_zero(itotv)=',result_zero(itotv),'result_first(itotv)=',result_first(itotv)

                endif


                if(nbackdT/=2.and.submodel<=0)then !20230216
                    if(iffix(itotv)==1) then ! displacement fixed

                        deltafi(itotv)=fixed(itotv)-result_zero(itotv)
                        !write(7,*)'itotv=',itotv,'fixed(itotv)=',fixed(itotv),'result_zero(itotv)=',result_zero(itotv)
                        if(ordert==1) then	 !x1
                            result_first(itotv)=(fixed(itotv)-result_zero(itotv))/ditime
                        else if(ordert==2) then
                            accmid=fixed(itotv)-(result_zero(itotv)+              &
                                ditime*result_first(itotv)+       &
                                .5*ditime**2*result_second(itotv))
                            accmid=accmid/(beeta2*ditime**2)
                            result_first(itotv)=result_first(itotv)+ditime*result_second(itotv) &
                                +beeta1*ditime*accmid
                            result_second(itotv)=result_second(itotv)+accmid
                        endif	 ! end x1
                        result_zero(itotv)=fixed(itotv)

                        !write(7,*)'itotv=',itotv,'result_zero(itotv)=',result_zero(itotv),'result_first(itotv)=',result_first(itotv)


                    else if(iffix(itotv)==2) then ! velocity fixed
                        if(ordert==1) then
                            deltafi(itotv)=ditime*(result_first(itotv)*(1-theta1)+theta1*fixed(itotv))
                        else if(ordert==2) then
                            accmid=fixed(itotv)-result_first(itotv)-ditime*result_second(itotv)
                            accmid=accmid/(beeta1*ditime)
                            deltafi(itotv)=ditime*result_first(itotv)+       &
                                ditime**2*(.5*result_second(itotv)+beeta2*accmid)
                            result_second(itotv)=result_second(itotv)+accmid
                        endif
                        result_zero(itotv)=result_zero(itotv)+deltafi(itotv)
                        result_first(itotv)=fixed(itotv)
                    else	if(iffix(itotv)==3) then					  !acceleration fixed
                        accmid=fixed(itotv)-result_second(itotv)
                        deltafi(itotv)=ditime*result_first(itotv)+       &
                            ditime**2*(.5*result_second(itotv)+beeta2*accmid)
                        result_zero(itotv)=result_zero(itotv)+deltafi(itotv)
                        result_first(itotv)=result_first(itotv)+ditime*result_second(itotv) &
                            +beeta1*ditime*accmid
                        result_second(itotv)=fixed(itotv)
                    endif
                endif  !20210324
            endif		 !end y1
        end do

    enddo

    delitfi=deltafi

    if(type_problem=='F'.and.ngapb/=0)then !!20121001
        do igapb=1,ngapb
            if(gapb(igapb)%nrdof==0)cycle
            gapb(igapb)%rdisp_zero=gapb(igapb)%rdisp_zero+ditime*gapb(igapb)%rdisp_first+       &
                .5*ditime**2*gapb(igapb)%rdisp_second
            gapb(igapb)%rdisp_delitfi=ditime*gapb(igapb)%rdisp_first+       &
                .5*ditime**2*gapb(igapb)%rdisp_second
            gapb(igapb)%rdisp_deltafi=gapb(igapb)%rdisp_delitfi
            gapb(igapb)%rdisp_first=gapb(igapb)%rdisp_first+.5*ditime*gapb(igapb)%rdisp_second
            !
            !
            !allocate(disl(gapb(igapb)%nrdof)) !20171130
            !disl=ditime*gapb(igapb)%rdisp_first+.5*ditime**2*gapb(igapb)%rdisp_second !此处为刚体转动加速度增量
            !   call dis_modify(igapb,disl)
            !disl=.5*ditime*gapb(igapb)%rdisp_second
            !   call vel_modify(igapb,disl)
            !   deallocate(disl) !20171130
            !

        end do   !igapb
    endif  !!20121001




    !! Simo rifai
    DO igroup =1,ngroup
        field1= group(igroup)%fieldid(1:1)
        special= group(igroup)%special
        index = group(igroup)%index
        if (appear(igroup)>0.and.field1=='U')then
            if ((index==5.or.index==9.or.index==16.or.index==18).and.special(1:1)=='B')then
                DO ielgroup = 1,group(igroup)%nelgroup
                    ielem = group(igroup)%list(ielgroup)
                    if (order_time_mdofn(1)==1)  then
                        element(ielem)%alfa=ditime*element(ielem)%alfa_first
                    elseif(order_time_mdofn(1)==2)  then
                        element(ielem)%alfa=ditime*element(ielem)%alfa_first
                        element(ielem)%alfa=ditime*element(ielem)%alfa_first+    &
                            .5*ditime*ditime*element(ielem)%alfa_second
                    endif
                end do
            endif
        endif
    end do
    !!simo rifai
    END  SUBROUTINE PREDICT

    SUBROUTINE varupdate
    integer(ink) idofn,ordert,ipoin,itotv,igapb,idimn,jdimn,i1,i2,ipoin1,ipoin2,igaps,ipairs,npairs
    integer(ink) igroup,ielem,ielgroup,index,ic,jdofn  !simo_rifai
    real(irk) xxac
    real(irk),allocatable::disl(:)
    character(5 )fieldid  !simo_rifai
    character(10)special  !simo_rifai

    ic=0
    if(type_load=='LOAD2'.and.idiv==2)ic=1
    if(type_load/='LOAD2')ic=1
    !write(7,*)'ipoin,idofn,itotv,result_zero(itotv)'
    !write(7,*)'iffix(nodfn(1:3,25))=',iffix(nodfn(1:3,25))
    !write(7,*)'result_zero(iffix(nodfn(1:3,25)))1=',result_zero(iffix(nodfn(1:3,25)))
    delitfi=0.
    if(block_stab>=1.and.ebody/=1) goto 10
    !allocate(disl(ndimn))
    do idofn=1,cdofn
        jdofn=lcdofn(idofn)               !20200226
        if(jdofn==10.and.outintr/=0)cycle !20200226
        ordert=order_time_mdofn(lcdofn(idofn))
        do ipoin=1,npoin
            itotv=nodfn(idofn,ipoin)
            if(itotv==0)cycle
            if (itotv/=0.and.(iffix(itotv)==0.or.trans(itotv)%nintf/=0))then
                !if (itotv/=0.and.iffix(itotv)==0)then
                if (ordert==0) then      !! for ordert


                    if(ic==1)&
                        result_zero(itotv)=result_zero(itotv)+result(itotv)
                    delitfi(itotv)    =result(itotv)
                    !write(7,*)'ipoin=',ipoin,'idofn=',idofn,'result=',result(itotv)
                elseif(ordert==1) then
                    if(ic==1)then
                        result_zero(itotv)=result_zero(itotv)+theta1*ditime*result(itotv)
                        result_first(itotv)=result_first(itotv)+result(itotv)
                        if(type_problem=='F'.and.lcdofn(idofn)==8)then   !20220721
                            result_second(itotv)=result_second(itotv)+result(itotv)/ditime
                        endif
                    endif
                    delitfi(itotv)    =theta1*ditime*result(itotv)
                    !write(7,20)ipoin,idofn,itotv,result_zero(itotv)

                elseif(ordert==2)then
                    if(ic==1)then
                        result_zero(itotv)=result_zero(itotv)+beeta2*ditime**2*result(itotv)
                        result_first(itotv)=result_first(itotv)+beeta1*ditime*result(itotv)
                        result_second(itotv)=result_second(itotv)+result(itotv)
                    endif
                    delitfi(itotv)=beeta2*ditime**2*result(itotv)
                    !write(7,20)ipoin,idofn,itotv,result_zero(itotv)


                endif !! for ordert

                ! if(idofn<=ndimn.or.idofn==7)deltafi(itotv)=deltafi(itotv)+delitfi(itotv)
                if(idofn<=7)&   !20220607
                    deltafi(itotv)=deltafi(itotv)+delitfi(itotv)
                !write(7,*)'itotv=',itotv,'deltafi0=',deltafi(itotv),'delitfi(itotv)=',delitfi(itotv)

            endif
        enddo     !! for ipoin
    end do     !! for idofn
20  format(3I10,e15.6)

    !deallocate(disl)
10  continue
    !write(7,*)'ngapb=',ngapb
    !write(7,*)'result_zero(iffix(nodfn(1:3,25)))2=',result_zero(iffix(nodfn(1:3,25)))

    if(ngapb/=0)then
        do igapb=1,ngapb
            if(gapb(igapb)%nrdof==0)cycle   !20231026
            allocate(disl(gapb(igapb)%nrdof))

            disl=gapb(igapb)%rdisp_inc !此处为刚体转动加速度增量
            write(7,*)'igapb=',igapb,'disl=',disl
            if(type_problem=='Q')then !!20121001
                gapb(igapb)%rdisp_delitfi=disl
                if(ic==1) &
                    gapb(igapb)%rdisp_zero=gapb(igapb)%rdisp_zero+disl
                gapb(igapb)%rdisp_deltafi=gapb(igapb)%rdisp_deltafi+disl
                !write(7,*)'before_dis_modify*****'
                call dis_modify(igapb,disl)
                write(7,*)'disl=',disl,'rdisp_zero=',gapb(igapb)%rdisp_zero
                write(7,*)'ext_force=',gapb(igapb)%ext_force
                write(7,*)'constrained stiffness='
                do idimn=1,3*(ndimn-1)
                    if(abs(gapb(igapb)%rdisp_zero(idimn))>1.e-20) &
                        write(7,30)idimn,gapb(igapb)%ext_force(idimn)/gapb(igapb)%rdisp_zero(idimn)
                end do
30              format(i10,e15.6)

            elseif(type_problem=='F')then !!20121001

                !write(7,*)'igapb=',igapb,'ic=',ic
                if(ic==1) &
                    gapb(igapb)%rdisp_zero=gapb(igapb)%rdisp_zero+beeta2*ditime**2*disl
                gapb(igapb)%rdisp_delitfi=beeta2*ditime**2*disl
                gapb(igapb)%rdisp_deltafi=gapb(igapb)%rdisp_deltafi+beeta2*ditime**2*disl
                call dis_modify(igapb,beeta2*ditime**2*disl)
                if(ic==1)then
                    gapb(igapb)%rdisp_first=gapb(igapb)%rdisp_first+beeta1*ditime*disl
                    call vel_modify(igapb,beeta1*ditime*disl)
                    gapb(igapb)%rdisp_second=gapb(igapb)%rdisp_second+disl
                    call acc_modify(igapb,disl)
                endif
                !call acc_check(igapb, gapb(igapb)%rdisp_second,acccheck)
            endif  !!20121001
            deallocate(disl)
        end do   !igapb
    end if

    !write(7,*)'result_zero(iffix(nodfn(1:3,25)))3=',result_zero(iffix(nodfn(1:3,25)))

    !deallocate(acccheck)

    !! Simo rifai
    if(ic==1)then
        DO igroup =1,ngroup
            fieldid= group(igroup)%fieldid
            special= group(igroup)%special
            index = group(igroup)%index
            if (appear(igroup)>0.and.fieldid(1:1)=='U')then
                if ((index==5.or.index==9.or.index==16.or.index==18).and.special(1:1)=='B')then
                    DO ielgroup = 1,group(igroup)%nelgroup
                        ielem = group(igroup)%list(ielgroup)
                        call updalfa(ielem,fieldid)
                        if (order_time_mdofn(1)==0)  then
                            element(ielem)%alfa=element(ielem)%alfa+element(ielem)%alfa_it
                        elseif(order_time_mdofn(1)==1)  then
                            element(ielem)%alfa=element(ielem)%alfa+theta1*ditime*element(ielem)%alfa_it
                            element(ielem)%alfa_first=element(ielem)%alfa_first+element(ielem)%alfa_it
                        elseif(order_time_mdofn(1)==2)  then
                            element(ielem)%alfa=element(ielem)%alfa+beeta2*ditime*ditime*element(ielem)%alfa_it
                            element(ielem)%alfa_first=element(ielem)%alfa_first+beeta1*ditime*element(ielem)%alfa_it
                            element(ielem)%alfa_second=element(ielem)%alfa_second+element(ielem)%alfa_it
                        endif
                    end do
                endif
            endif
        end do
    endif

    END SUBROUTINE varupdate

    SUBROUTINE relative_dis_watertight   !20231007 止水
    character (10) model,field1,material
    integer(ink) ielem,nnode,nevab,ngaus,nnode_half,inode,idofn,ielgroup,igroup,idimn,index
    integer(ink),allocatable:: lnods(:),ldofs(:)
    real   (irk),allocatable:: eldis(:),rot(:,:),shapecg(:,:),relat_dis_nod(:,:),nordis(:,:),relat_dis_gaus(:,:)

    DO igroup =1,ngroup
        field1= group(igroup)%fieldid(1:1)
        if(appear(igroup)>0.and.field1=='U')then
            matno = group(igroup)%matno
            material=props(matno)%mechanical%solid%material

            index    =group(igroup)%index
            nnode    =elkn(index)%el_field(1)%nnode_f

            if(material/='GOODMAN') cycle
            model=props(matno)%mechanical%solid%Goodman%model
            if(model/='WATERTIGHT') cycle
            jndex=1
            if(ndimn==3)jndex=5
            order_int=elkn(jndex)%el_field(1)%order_intrules(1)
            ngaus=elkn(jndex)%ggaus(order_int)%ngaus
            nevab=nnode*ndimn

            allocate(lnods(nnode),ldofs(nevab),eldis(nevab),nordis(ndimn,nnode),rot(ndimn,ndimn))
            nnode_half=nnode/2
            allocate(shapecg(nnode_half,ngaus),relat_dis_gaus(ndimn,ngaus),relat_dis_nod(ndimn,nnode))

            DO ielgroup = 1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                lnods = element(ielem)%field(1)%lnods_f
                ldofs = element(ielem)%field(1)%ldofs_f
                eldis = deltafi(ldofs)

                rot=element(ielem)%rotation

                do inode=1,nnode
                    idofn=(inode-1)*ndimn
                    nordis(:,inode)=rot.x.eldis(idofn+1:idofn+ndimn)
                end do

                if(ndimn==3) then
                    do inode=1,nnode_half
                        relat_dis_nod(:,inode)=nordis(:,inode+nnode_half)-nordis(:,inode)
                        relat_dis_nod(:,inode+nnode_half)=relat_dis_nod(:,inode)
                    end do
                else  ! 2D
                    relat_dis_nod(:,1)=nordis(:,4)-nordis(:,1)
                    relat_dis_nod(:,4)=relat_dis_nod(:,1)
                    relat_dis_nod(:,2)=nordis(:,3)-nordis(:,2)
                    relat_dis_nod(:,3)=relat_dis_nod(:,2)
                endif
                element(ielem)%field(1)%relat_dis_nod=element(ielem)%field(1)%relat_dis_nod0  &
                    +relat_dis_nod

                shapecg = elkn(jndex)%ggaus(order_int)%shape(:,:)
                do idimn=1,ndimn
                    relat_dis_gaus(idimn,:)=transpose(shapecg).x.relat_dis_nod(idimn,:)
                end do

                element(ielem)%field(1)%relat_dis_gaus=element(ielem)%field(1)%relat_dis_gaus0  &
                    +relat_dis_gaus

            end do !ielgroup

            deallocate(lnods,ldofs,eldis,nordis,rot,shapecg,relat_dis_nod,relat_dis_gaus)
        endif
    end do !igroup

    END SUBROUTINE relative_dis_watertight

    SUBROUTINE construction_dis_modify
    character(1)field1
    integer(ink) igroup,ielgroup,ielem,idofn
    integer(ink),pointer::lnods(:)
    integer(ink),allocatable::icd(:)
    real   (irk) zz,hh,coef

    allocate(icd(npoin))
    icd=0

    DO igroup =1,ngroup
        field1= group(igroup)%fieldid(1:1)
        if(field1/='U')cycle
        if (appear_process(igroup,iblks-1)==0.and.appear_process(igroup,iblks)/=0)then
            do ielgroup =1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                lnods=>element(ielem)%field(1)%lnods_f
                icd(lnods)=1
                nullify(lnods)
            end do
        end if
    end do


    do igroup =1,ngroup
        field1= group(igroup)%fieldid(1:1)
        if(field1/='U')cycle
        if (appear_process(igroup,iblks-1)/=0.and.appear_process(igroup,iblks)/=0)then
            DO ielgroup =1,group(igroup)%nelgroup
                ielem = group(igroup)%list(ielgroup)
                lnods=>element(ielem)%field(1)%lnods_f
                icd(lnods)=0
                nullify(lnods)
            end do
        end if
    end do

    if(iblks>1) then
        hh= hdam(iblks)-hdam(iblks-1)
    else
        hh= hdam(iblks)
    endif

    !if(iblks==12)write(7,*)'hh=',hh
    do ipoin=1,npoin
        if(icd(ipoin)==0)cycle
        zz=hdam(iblks)-coord(ndimn,ipoin)
        if(zz<0.)zz=0.
        coef=2.*zz/(hh+zz)
        !if(iblks==12)then
        !   write(7,*)'ipoin=',ipoin,'coef=',coef
        !endif
        do idofn=1,cdofn
            itotv=nodfn(idofn,ipoin)
            if(itotv==0)cycle
            result_zero(itotv)=result_zero(itotv)*coef
        enddo
    end do

    deallocate(icd)
    end SUBROUTINE construction_dis_modify

    SUBROUTINE acc_modify(igapb,disl1)
    integer(ink) igapb,ipoin,jpoin,idimn,jdimn,itotv,kdimn
    real   (irk) xxac,disl1(:)
    if(gapb(igapb)%nrdof==0)return
    kdimn=ndimn
    if(block_stab==1)kdimn=3*(ndimn-1)
    do jpoin=1,gapb(igapb)%npblock
        ipoin=gapb(igapb)%nodeblock(jpoin)
        do idimn=1,kdimn  !idimn
            itotv=nodfn(idimn,ipoin)
            if(itotv==0)cycle
            xxac=0.
            do jdimn=1,gapb(igapb)%nrdof
                xxac=xxac+disl1(jdimn)*gapb(igapb)%npdisp(idimn,jpoin,jdimn)
            end do
            result_second(itotv)=result_second(itotv)+xxac
        end do   !idimn
    end do   !jpoin
    end SUBROUTINE acc_modify

    SUBROUTINE acc_check(igapb,disl1,acccheck)
    integer(ink) igapb,ipoin,jpoin,idimn,jdimn,itotv
    real   (irk) xxac,disl1(:),acccheck(:)

    if(gapb(igapb)%nrdof==0)return
    do jpoin=1,gapb(igapb)%npblock
        ipoin=gapb(igapb)%nodeblock(jpoin)
        do idimn=1,ndimn  !idimn
            itotv=nodfn(idimn,ipoin)
            if(itotv==0)cycle
            xxac=0.
            do jdimn=1,gapb(igapb)%nrdof
                xxac=xxac+disl1(jdimn)*gapb(igapb)%npdisp(idimn,jpoin,jdimn)
            end do
            acccheck(itotv)=acccheck(itotv)+xxac
        end do   !idimn
    end do   !jpoin
    end SUBROUTINE acc_check

    SUBROUTINE acc_rigid
    integer(ink) igapb,ipoin,jpoin,idimn,jdimn,itotv,kdimn
    real   (irk) xxac
    real   (irk), allocatable::disl(:)
    !write(7,*)'igapb=',igapb
    kdimn=ndimn
    if(block_stab==1)kdimn=3*(ndimn-1)
    do igapb=1,ngapb
        if(gapb(igapb)%nrdof==0)cycle

        allocate(disl(gapb(igapb)%nrdof))
        disl=gapb(igapb)%rdisp_second
        !   disl=100.
        do jpoin=1,gapb(igapb)%npblock
            ipoin=gapb(igapb)%nodeblock(jpoin)
            do idimn=1,kdimn  !idimn
                itotv=nodfn(idimn,ipoin)
                if(itotv==0)cycle
                xxac=0.
                do jdimn=1,gapb(igapb)%nrdof
                    xxac=xxac+disl(jdimn)*gapb(igapb)%npdisp(idimn,jpoin,jdimn)
                end do
                !   write(7,*)'itotv=',itotv,'acc=',xxac
                result_second(itotv)=xxac
            end do   !idimn
        end do   !jpoin
        deallocate(disl)
    end do
    end SUBROUTINE acc_rigid

    SUBROUTINE vel_modify(igapb,disl1)
    integer(ink) igapb,ipoin,jpoin,idimn,jdimn,itotv,kdimn
    real   (irk) xxac,disl1(:)

    kdimn=ndimn
    if(block_stab==1)kdimn=3*(ndimn-1)
    if(gapb(igapb)%nrdof==0)return
    do jpoin=1,gapb(igapb)%npblock
        ipoin=gapb(igapb)%nodeblock(jpoin)
        do idimn=1,kdimn  !idimn
            itotv=nodfn(idimn,ipoin)
            if(itotv==0)cycle
            xxac=0.
            do jdimn=1,gapb(igapb)%nrdof
                xxac=xxac+disl1(jdimn)*gapb(igapb)%npdisp(idimn,jpoin,jdimn)
            end do
            result_first(itotv)=result_first(itotv)+xxac
        end do   !idimn
    end do   !jpoin
    end SUBROUTINE vel_modify

    SUBROUTINE dis_modify(igapb,disl1)
    integer(ink) igapb,ipoin,jpoin,idimn,jdimn,itotv,kdimn
    real   (irk) xxac,disl1(:)

    kdimn=ndimn
    if(block_stab==1)kdimn=3*(ndimn-1)
    if(gapb(igapb)%nrdof==0)return
    !write(7,*)'igapb=',igapb,'disl1=',disl1,'gapb(igapb)%npblock=',gapb(igapb)%npblock
    do jpoin=1,gapb(igapb)%npblock
        ipoin=gapb(igapb)%nodeblock(jpoin)
        do idimn=1,kdimn  !idimn
            itotv=nodfn(idimn,ipoin)
            if(itotv==0)cycle
            xxac=0.
            do jdimn=1,gapb(igapb)%nrdof
                xxac=xxac+disl1(jdimn)*gapb(igapb)%npdisp(idimn,jpoin,jdimn)
            end do
            !if (iffix(itotv)==0.or.iffix(itotv)==4)then
            delitfi(itotv)=delitfi(itotv)+xxac
            deltafi(itotv)=deltafi(itotv)+xxac
            result_zero(itotv)=result_zero(itotv)+xxac
            if(nbackf/=0)result_zero_g(itotv)=result_zero_g(itotv)+xxac
            !if(itotv>=100.and.itotv<=120) &
            !write(7,*)'itotv=',itotv,'delitfi=',delitfi(itotv),'deltafi=',deltafi(itotv),'result_zero=',result_zero(itotv)
            !else
            !write(7,*)'ipoin=',ipoin,'idimn=',idimn,'itotv=',itotv,'iffix=',iffix(itotv),'xxac=',xxac,'result_zero(itotv)=',result_zero(itotv)
            !endif

        end do   !idimn
    end do   !jpoin
    end SUBROUTINE dis_modify

    SUBROUTINE varupdate_w !freq2006

    integer(ink) idofn,ordert,ipoin,itotv

    !delitfi=0.
    do idofn=1,cdofn
        do ipoin=1,npoin
            itotv=nodfn(idofn,ipoin)
            if (itotv/=0.and.iffix(itotv)==0)then
                result_zero(itotv)=abs(resultw(itotv))*ttime**2
            elseif(itotv/=0.and.iffix(itotv)/=0)then
                result_zero(itotv)=fixed(itotv)*ttime**2
                resultw(itotv)=fixed(itotv)
            endif
            !        if (itotv/=0)then
            !           delitfi(itotv)=result_zero(itotv)
            !           if(idofn<=7)deltafi(itotv)=deltafi(itotv)+delitfi(itotv)
            !        endif
        enddo     !! for ipoin
    end do     !! for idofn

    END SUBROUTINE varupdate_w

    subroutine updalfa(ielem,fieldid)
    !
    character(5) fieldid
    integer(ink) ielem,nevab,aevab,np
    integer(ink), pointer::ldofsp(:),ldofs(:)
    real   (irk), pointer::rh(:),estift(:,:),estifhi(:,:), qmatxa(:,:)
    real   (irk), allocatable::eldis(:),rh0(:,:),alfa0(:,:),rh1(:),elpw(:)

    ldofs => element(ielem)%field(1)%ldofs_f
    nevab=size(ldofs)
    aevab=size(element(ielem)%rh)
    if (fieldid(1:2)=='UW') then
        ldofsp => element(ielem)%field(2)%ldofs_f
        np = size(ldofsp)
        allocate(elpw(np))
        elpw  = delitfi(ldofsp)
    endif
    allocate(eldis(nevab))
    rh     =>element(ielem)%rh
    eldis = delitfi(ldofs )

    estift =>element(ielem)%estift
    estifhi=>element(ielem)%estifh
    allocate(rh0(aevab,1),alfa0(aevab,1),rh1(aevab))
    rh1=estift.x.eldis

    if(fieldid(1:2)=='UW') then
        qmatxa => element(ielem)%qmatxa
        rh1 = rh1 - ( qmatxa.x.elpw )
    endif

    rh0(:,1)=-rh-rh1

    call householder(estifhi,rh0,alfa0) !3
    !
    element(ielem)%alfa_it=alfa0(:,1)

    nullify(rh,ldofs,estift,estifhi)
    deallocate(alfa0,rh0,rh1,eldis)

    if (fieldid(1:2)=='UW') then
        nullify(ldofsp)
        deallocate(elpw)
    endif

    end  subroutine updalfa

    SUBROUTINE ALGORT
    !
    integer(ink) KOUNT

    KRESL=0
    IF(type_nl.EQ.1.and.istep==inc_step.AND.IITER.EQ.1)   KRESL=1
    IF(type_nl.EQ.2.AND.IITER.le.2)                       KRESL=1
    IF(type_nl.EQ.4)                                      KRESL=1    ! Full newton
    IF(type_nl.EQ.5.AND.IITER.EQ.1)                       KRESL=1
    IF(type_nl.EQ.6.AND.IITER.EQ.2)                       KRESL=1
    IF(type_nl.EQ.7.AND.IITER.ge.2)                       KRESL=1
    if(type_nl.eq.8.and.iiter==1)                         kresl=1
    IF(type_nl==9.AND.iblks==(lblks+1).and. &
        iincs==1.and.istep==inc_step.and.IITER.EQ.1)          kresl=1
    if(type_nl.eq.10)                                     kresl=0
    !   if(kstat==2.and.iiter<=2)                             kresl=1
    if(gamamax/=0)then !20231215YL
        if(istep/=1)kresl=0
    endif !20231215YL



    if (nlayer==2) then
        kresl_layer1=0
        IF(type_nl_layer1.EQ.1.AND.istep==inc_step.AND.IITER.EQ.1)       kresl_layer1=1
        IF(type_nl_layer1.EQ.2.AND.IITER.le.2)                           kresl_layer1=1
        IF(type_nl_layer1.EQ.4)                                          kresl_layer1=1  ! Full newton
        IF(type_nl_layer1.EQ.5.AND.IITER.EQ.1)                           kresl_layer1=1
        IF(type_nl_layer1==9.AND.  &
            iblks==(lblks+1).and.iincs==1.and.istep==inc_step.and.IITER.EQ.1)kresl_layer1=1
        kresl_layer2=0
        IF(type_nl_layer2.EQ.1.AND.istep==inc_step.AND.IITER.EQ.1)       kresl_layer2=1
        IF(type_nl_layer2.EQ.2.AND.IITER.le.2)                           kresl_layer2=1
        IF(type_nl_layer2.EQ.4)                                          kresl_layer2=1  ! Full newton
        IF(type_nl_layer2.EQ.5.AND.IITER.EQ.1)                           kresl_layer2=1
        IF(type_nl_layer2==9.AND.  &
            iblks==(lblks+1).and.iincs==1.and.istep==inc_step.and.IITER.EQ.1)kresl_layer2=1
        kresl=0
        if(kresl_layer1/=0.or.kresl_layer2/=0)                           kresl=1
    endif
    !
    KMASS=0
    if (type_problem/='Q'.or.stabpw==1) then
        IF (NMASS.EQ.0.OR.(istep==inc_step.AND.IITER.EQ.1)) THEN
            KMASS=1
        ELSE
            IF (IITER.EQ.1) THEN
                KOUNT=(istep/NMASS)*NMASS
                IF(KOUNT.EQ.istep)                                          KMASS=1
            END IF
        END IF
    endif
    !
    KSMAT=0
    IF (NSMAT.EQ.0.OR.(istep==inc_step.AND.IITER.EQ.1)) THEN
        KSMAT=1
    ELSE
        IF (IITER.EQ.1) THEN
            KOUNT=(istep/NSMAT)*NSMAT
            IF(KOUNT.EQ.istep)                                             KSMAT=1
        END IF
    END IF
    if(type_problem=='Q'.or.outinp>0)                                                        ksmat=0  !20220623
    !
    KHMAT=0
    IF (NHMAT.EQ.0.OR.(istep==inc_step.AND.IITER.EQ.1)) THEN
        KHMAT=1
    ELSE
        IF (IITER.EQ.1) THEN
            KOUNT=(istep/NHMAT)*NHMAT
            IF(KOUNT.EQ.istep)                                             KHMAT=1
        END IF
    END IF
    if(outinp>0)                                                        khmat=0  !20220623

    !   temperature

    KTSMAT=0
    IF (NTSMAT.EQ.0.OR.(istep==inc_step.AND.IITER.EQ.1)) THEN
        KTSMAT=1
    ELSE
        IF (IITER.EQ.1) THEN
            KOUNT=(istep/NTSMAT)*NTSMAT
            IF(KOUNT.EQ.istep)                                             KTSMAT=1
        END IF
    END IF
    !if(outintr.ne.0.or.outinp/=0)                                       ktsmat=0  !20220623
    if(outintr.ne.0)                                       ktsmat=0   !20220623
    !
    KTHMAT=0
    IF (NTHMAT.EQ.0.OR.(istep==inc_step.AND.IITER.EQ.1)) THEN
        KTHMAT=1
    ELSE
        IF (IITER.EQ.1) THEN
            KOUNT=(istep/NTHMAT)*NTHMAT
            IF(KOUNT.EQ.istep)                                             KTHMAT=1
        END IF
    END IF
    !if(outintr.ne.0.or.outinp/=0)                                       kthmat=0 !20220623
    if(outintr.ne.0)                                       kthmat=0 !20220623


    !   temperature
    !
    KQMAT=0
    IF (NQMAT.EQ.0.OR.(istep==inc_step.AND.IITER.EQ.1)) THEN
        KQMAT=1
    ELSE
        IF (IITER.EQ.1) THEN
            KOUNT=(istep/NQMAT)*NQMAT
            IF(KOUNT.EQ.istep)                                             KQMAT=1
        END IF
    END IF
    !
    KSWKW=0
    IF (NSWKW.EQ.0.or.(istep.eq.inc_step.and.iiter.eq.1)) THEN
        KSWKW=1
    ELSE
        IF (IITER.EQ.1) THEN
            KOUNT=(ISTEP/NSWKW)*NSWKW
            IF(KOUNT.EQ.ISTEP)                                             KSWKW=1
        END IF
    END IF

    !
    kldfl=0
    IF (NLDFL.EQ.0.OR.(istep==inc_step.AND.IITER.EQ.1)) THEN
        KLDFL=1
    ELSE
        IF (IITER.EQ.1) THEN
            KOUNT=(istep/NLDFL)*NLDFL
            IF(KOUNT.EQ.istep)                                             KLDFL=1
        END IF
    END IF
    !
    kgrav=0
    IF (NGRAV.EQ.0.OR.(istep==inc_step.AND.IITER.EQ.1)) THEN
        KGRAV=1
    ELSE
        IF (IITER.EQ.1) THEN
            KOUNT=(istep/NGRAV)*NGRAV
            IF(KOUNT.EQ.istep)                                             KGRAV=1
        END IF
    END IF


    END SUBROUTINE ALGORT
