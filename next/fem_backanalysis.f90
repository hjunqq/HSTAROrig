    subroutine parameter_back_analysis(N,M) !20190810
    !N随机变量数,M测点数
    implicit none

    INTEGER  i,j, N, M,ivalue,inode,idofn,jnode,jtotv,iobstimes,  &
        iblks_bp,iincs_bp,istep_bp,iobse_bp,jvalue,ivalue_point
    !** SOLUTION VECTOR. CONTAINS VALUES X FOR F(X)
    DOUBLE PRECISION    X (N)
    DOUBLE PRECISION,allocatable:: errn(:,:),fjac2(:,:),errx(:,:),fjac1(:,:), &
        errm0(:,:),sigma(:),errn0(:,:),errx0(:,:)


    x=xvalue
    !call parameter_back_analysis_read
    IF (DTRNLSP_INIT (HANDLE,N,M,X,EPS, ITER1, ITER2, RS).NE. TR_SUCCESS) THEN
        !** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        PRINT *, '| ERROR IN DTRNLSP_INIT'
        !** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
        !** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
        !** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
        CALL MKL_FREE_BUFFERS
        !** AND STOP
        STOP 1
    END IF

    !   ** CHECKS THE CORRECTNESS OF HANDLE AND ARRAYS CONTAINING JACOBIAN MATRIX,
    !** OBJECTIVE FUNCTION, LOWER AND UPPER BOUNDS, AND STOPPING CRITERIA.
    IF (DTRNLSP_CHECK (HANDLE, N,M, FJAC, FVEC, EPS, INFO) .NE. TR_SUCCESS) THEN
        !** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        PRINT *, '| ERROR IN DTRNLSPBC_INIT'
        !** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
        !** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
        !** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
        CALL MKL_FREE_BUFFERS
        !** AND STOP
        STOP 1
    ELSE
        !write(7,*)'size(fjac)=',size(fjac),'size(fvec)=',size(fvec)
        !write(7,*)'info=',info(1:4)
        !** THE HANDLE IS NOT VALID.
        IF( INFO(1) .NE. 0 .OR.    &
            !** THE FJAC ARRAY IS NOT VALID.
            INFO(2) .NE. 0 .OR.   &
            !** THE FVEC ARRAY IS NOT VALID.
            INFO(3) .NE. 0 .OR. &
            !** THE EPS ARRAY IS NOT VALID.
            INFO(4) .NE. 0 ) THEN
            PRINT *, '| INPUT PARAMETERS ARE NOT VALID'
            !** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
            !** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
            !** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
            CALL MKL_FREE_BUFFERS
            !** AND STOP
            STOP 1
        END IF
    END IF

    !** SET INITIAL RCI CYCLE VARIABLES
    RCI_REQUEST = 0
    SUCCESSFUL = 0



    DO WHILE (SUCCESSFUL == 0)
        !** CALL TR SOLVER
        !**   HANDLE        IN/OUT: TR SOLVER HANDLE
        !**   FVEC          IN:     VECTOR
        !**   FJAC          IN:     JACOBI MATRIX
        !**   RCI_REQUEST   IN/OUT: RETURN NUMBER WHICH DENOTE NEXT STEP FOR PERFORMING
        IF (DTRNLSP_SOLVE (HANDLE, FVEC, FJAC, RCI_REQUEST)  &
            .NE. TR_SUCCESS) THEN
            !** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
            PRINT *, '| ERROR IN DTRNLSP_SOLVE'
            !** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
            !** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
            !** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
            CALL MKL_FREE_BUFFERS
            !** AND STOP
            STOP 1
        END IF
        !** ACCORDING WITH RCI_REQUEST VALUE WE DO NEXT STEP

        res=DTRNLSP_GET (HANDLE, ITER, ST_CR, R1, R2)
        write(7,*)'res=',res
        write(7,*)'r1=',r1,'r2=',r2
        write(7,*)'RCI_REQUEST=',RCI_REQUEST
        write(7,*)'iter=',iter
        write(7,*)'xvalue=',x


        !print *,'res=',res
        !print *,'r1=',r1,'r2=',r2
        !print *,'RCI_REQUEST=',RCI_REQUEST
        print *,'iter=',iter
        print *,'xvalue=',x



        SELECT CASE (RCI_REQUEST)
        CASE (-1, -2, -3, -4, -5, -6)
            !**   STOP RCI CYCLE
            SUCCESSFUL = 1
        CASE (1)
            !**   RECALCULATE FUNCTION VALUE
            !**     M               IN:     DIMENSION OF FUNCTION VALUE
            !**     N               IN:     NUMBER OF FUNCTION VARIABLES
            !**     X               IN:     SOLUTION VECTOR
            !**     FVEC            OUT:    FUNCTION VALUE F(X)
            xvalue=x
            do i=1,n
                if(para_back(i)%mode_transform==0)then
                    xvalue(i)=xvalue(i)*para_back(i)%factor
                elseif(para_back(i)%mode_transform==1)then
                    xvalue(i)=para_back(i)%factor/xvalue(i)
                endif
            end do
            ttime=0.
            trstep=0   !20230430
            Value_observ(:)%value_computation=0.

            call process_analysis  !20190810
            fvec=0.
            !write(7,*)'ivalue,jvalue,value_mesure,value_computation'
            do ivalue=1,m
                if(Value_observ(ivalue)%ic==0)cycle
                Fvec(ivalue)=Value_observ(ivalue)%value_measure-Value_observ(ivalue)%value_computation
                jvalue=Value_observ(ivalue)%jvalue
                if(jvalue/=0) &   !20230709
                    Fvec(ivalue)=Fvec(ivalue)+Value_observ(jvalue)%value_computation
                ivalue_point=Value_observ(ivalue)%ivalue_point

                !if(ivalue_point==169.or.ivalue_point==170)then
                idofn =lmdofn(Value_observ(ivalue)%idofn)
                !write(7,*)'ivalue_point=',ivalue_point,'idofn=',idofn
                ! if(jvalue==0) & !20230717
                !write(7,*)ivalue,jvalue,Value_observ(ivalue)%value_measure,Value_observ(ivalue)%value_computation
                ! if(jvalue/=0) & !20230717
                !write(7,*)ivalue,jvalue,Value_observ(ivalue)%value_measure,  &
                !    (Value_observ(ivalue)%value_computation-Value_observ(jvalue)%value_computation)
                !endif

            end do
            !write(7,*)'fvec=',fvec
        CASE (2)
            !**   COMPUTE JACOBI MATRIX
            !**     EXTENDED_POWELL IN:     EXTERNAL OBJECTIVE FUNCTION
            !**     N               IN:     NUMBER OF FUNCTION VARIABLES
            !**     M               IN:     DIMENSION OF FUNCTION VALUE
            !**     FJAC            OUT:    JACOBI MATRIX
            !**     X               IN:     SOLUTION VECTOR
            !**     JAC_EPS         IN:     JACOBI CALCULATION PRECISION

            if(balgor==0)then
                do  i=1, N
                    xvalue=x
                    Xvalue(i)=X(i)-para_back(i)%factor_inc*X(i)

                    do j=1,n
                        if(para_back(j)%mode_transform==0)then
                            xvalue(j)=xvalue(j)*para_back(j)%factor
                        elseif(para_back(j)%mode_transform==1)then
                            xvalue(j)=para_back(j)%factor/xvalue(j)
                        endif
                    end do
                    ttime=0.
                    trstep=0
                    Value_observ(:)%value_computation=0.
                    call process_analysis  !20190810
                    fvec1=0.
                    do ivalue=1,m
                        if(Value_observ(ivalue)%ic==0)cycle
                        Fvec1(ivalue)=Value_observ(ivalue)%value_computation
                    end do
                    xvalue=x
                    Xvalue(i)=X(i)+para_back(i)%factor_inc*X(i)
                    do j=1,n
                        if(para_back(j)%mode_transform==0)then
                            xvalue(j)=xvalue(j)*para_back(j)%factor
                        elseif(para_back(j)%mode_transform==1)then
                            xvalue(j)=para_back(j)%factor/xvalue(j)
                        endif
                    end do
                    ttime=0.
                    trstep=0
                    Value_observ(:)%value_computation=0.
                    call process_analysis  !20190810
                    fvec2=0.
                    !write(7,*)'ivalue,Value_observ(ivalue)%ic,Fvec1(ivalue),Fvec2(ivalue)='
                    do ivalue=1,m
                        if(Value_observ(ivalue)%ic==0)cycle
                        Fvec2(ivalue)=Value_observ(ivalue)%value_computation
                        !write(7,15)ivalue,Value_observ(ivalue)%ic,Fvec1(ivalue),Fvec2(ivalue)
                    end do
                    FJAC(:,i)=-(Fvec2-Fvec1)/(2.*para_back(i)%factor_inc*x(i))
                    !15 format(2i5,2e15.5)
                    !     write(7,*)'i=',i
                    !write(7,*)'fjac=',fjac(:,i)
                end do
                !stop
            else if(balgor==1)then  !20220108
                do  i=1, N
                    do J=1, mvalue
                        FJAC(j,i)=-Value_observ(j)%dudx(i)
                    enddo
                enddo
            endif

            !                IF (DJACOBI (measure_computation, N, M, FJAC, X, JAC_EPS)   &
            !                     .NE. TR_SUCCESS) THEN
            !!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
            !                    PRINT *, '| ERROR IN DJACOBI'
            !!** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
            !!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
            !!** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
            !                    CALL MKL_FREE_BUFFERS
            !!** AND STOP
            !                    STOP 1
            !                     END IF
        ENDSELECT
    END DO

    allocate(fjac1(npoints_pbx,n),errn(n,mobstimes),fjac2(n,n),errx(n,mobstimes))

    allocate(errm0(npoints_pbx,mobstimes),sigma(n))

    errm0=0.
    !write(7,*)'m=',m
    do ivalue=1,m
        if(Value_observ(ivalue)%ic==0)cycle
        iblks_bp=Value_observ(ivalue)%iblks
        iincs_bp=Value_observ(ivalue)%iincs
        istep_bp=Value_observ(ivalue)%istep
        iobse_bp=Value_observ(ivalue)%iobse
        !print *,'ivalue=',ivalue,'iblks_bp=',iblks_bp,'iincs_bp=',iincs_bp,'istep_bp=',istep_bp


        jtotv=Value_observ(ivalue)%ivalue_point
        !print *,'jtotv=',jtotv,'iobse_bp=',iobse_bp

        iobstimes=para_block(iblks_bp)%para_nincs(iincs_bp)%tstep_bp(istep_bp,iobse_bp)
        errm0(jtotv,iobstimes)=fvec(ivalue)
    end do
    !write(7,*)'errm0='
    !      do i=1,mobstimes
    !     write(7,10)errm0(:,i)
    !     end do

10  format(10e15.5)

    do i=1,mobstimes
        fjac1=0.
        do ivalue=1,m
            if(Value_observ(ivalue)%ic==0)cycle
            jtotv=Value_observ(ivalue)%ivalue_point
            iblks_bp=Value_observ(ivalue)%iblks
            iincs_bp=Value_observ(ivalue)%iincs
            istep_bp=Value_observ(ivalue)%istep
            iobse_bp=Value_observ(ivalue)%iobse
            iobstimes=para_block(iblks_bp)%para_nincs(iincs_bp)%tstep_bp(istep_bp,iobse_bp)
            if(i==iobstimes)then
                !write(7,*),'i=',i,'ivalue=',ivalue
                fjac1(jtotv,:)=fjac(ivalue,:)
            endif
        end do
        !     write(7,*)'iobstimes=',i,'fjac1='
        !     do jtotv=1,npoints_pbx
        !     write(7,10)fjac1(jtotv,:)
        !     end do
        !
        !stop

        do jtotv=1,npoints_pbx
            errn(:,i)=transpose(fjac1).x.errm0(:,i)
        end do
        fjac2=transpose(fjac1).x.fjac1
        allocate(errn0(n,1),errx0(n,1))
        errn0(:,1)=errn(:,i)
        call householderx(fjac2,errn0,errx0)
        errn(:,i)=errn0(:,1)
        errx(:,i)=errx0(:,1)
        deallocate(errn0,errx0)
    end do

    !write(7,*)'errx='
    !    do i=1,mobstimes
    !   write(7,10)errx(:,i)
    !   end do
    !write(7,*)'errn='
    !    do i=1,mobstimes
    !   write(7,10)errn(:,i)
    !   end do
    !

    do j=1,n
        sigma(j)=0.

        do i=1,mobstimes
            sigma(j)=sigma(j)+errx(j,i)**2
        end do
        sigma(j)=sigma(j)/mobstimes
        sigma(j)=sqrt(sigma(j))
        !sigma(j)=sigma(j)/x(j)
    end do


    write(7,*)'sigma=',sigma
    write(7,*)'mobstimes=',mobstimes
    write(7,*)'x=',x

    deallocate(errn,fjac2,fjac1,errx,sigma,errm0)


    !** GET SOLUTION STATUSES
    !**   HANDLE            IN: TR SOLVER HANDLE
    !**   ITER              OUT: NUMBER OF ITERATIONS
    !**   ST_CR             OUT: NUMBER OF STOP CRITERION
    !**   R1                OUT: INITIAL RESIDUALS
    !**   R2                OUT: FINAL RESIDUALS
    IF (DTRNLSP_GET (HANDLE, ITER, ST_CR, R1, R2)  &
        .NE. TR_SUCCESS) THEN
        !** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        PRINT *, '| ERROR IN DTRNLSP_GET'
        !** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
        !** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
        !** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
        CALL MKL_FREE_BUFFERS
        !** AND STOP
        STOP 1
    END IF
    !** FREE HANDLE MEMORY
    IF (DTRNLSP_DELETE (HANDLE) .NE. TR_SUCCESS) THEN
        !** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        PRINT *, '| ERROR IN DTRNLSP_DELETE'
        !** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
        !** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
        !** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
        CALL MKL_FREE_BUFFERS
        !** AND STOP
        STOP 1
    END IF

    !** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
    !** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
    !** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
    CALL MKL_FREE_BUFFERS
    !** IF FINAL RESIDUAL LESS THEN REQUIRED PRECISION THEN PRINT PASS
    !        IF (R2 .LT. 1.D-5) THEN
    !            PRINT *, '|         DTRNLSP POWELL............PASS'
    !            STOP 0
    !!** ELSE PRINT FAILED
    !        ELSE
    !            PRINT *, '|         DTRNLSP POWELL............FAILED'
    !            STOP 1
    !        END IF


    end  subroutine parameter_back_analysis

    subroutine trust_region_back_analysis(N,M) !20220108
    implicit none

    integer(ink)  i,j,k, N, M,ivalue,inode,idofn,jnode,jtotv,iobstimes,  &
        iblks_bp,iincs_bp,istep_bp,iobse_bp,iter_tr,mtter,jvalue
    !** SOLUTION VECTOR. CONTAINS VALUES X FOR F(X)
    real(irk)   X(N),x0(N),x1,x2,alfak,deltak, &
        eta1,eta2,gama1,gama2,eps,delta0,eta01,eta02,rk, &
        deltab,dkstar,dkpre
    real(irk),allocatable:: errn(:,:),fjac2(:,:),errx(:,:),fjac1(:,:), &
        errm0(:,:),sigma(:),errn0(:,:),errx0(:,:),gk0(:),bk0(:,:)
    real(irk),allocatable::gk(:),gk1(:,:),Bk(:,:),xnew(:),funx(:),  &
        Qs(:),sk(:),yk(:),funx0(:),bs(:)

    allocate(gk(N),gk1(N,1),Bk(N,N),xnew(N),funx0(n),bs(n))
    eta1=trustp(1)%eta1;   eta2=trustp(1)%eta2
    gama1=trustp(1)%gama1; gama2=trustp(1)%gama2
    eps=trustp(1)%eps;    deltak=trustp(1)%delta0
    eta01=trustp(1)%eta01; eta02=trustp(1)%eta02
    mtter=trustp(1)%mtter; deltab=trustp(1)%deltab
    allocate(funx(mtter),Qs(mtter),gk0(N),sk(N),yk(N),bk0(N,N))

    x0=xvalue
    xnew=xvalue
    Bk=0.
    do i=1,n
        bk(i,i)=1.
    end do

    iter_tr=1

    x=xnew
    do i=1,n
        if(para_back(i)%mode_transform==0)then
            xvalue(i)=x(i)*para_back(i)%factor
        elseif(para_back(i)%mode_transform==1)then
            xvalue(i)=para_back(i)%factor/x(i)
        endif
    end do
    ttime=0.
    Value_observ(:)%value_computation=0.

    call process_analysis
    fvec=0.

    write(7,*)'ivalue,value_mesure,value_computation'
    do ivalue=1,m
        if(Value_observ(ivalue)%ic==0)cycle
        Fvec(ivalue)=Value_observ(ivalue)%value_measure-Value_observ(ivalue)%value_computation
        jvalue=Value_observ(ivalue)%jvalue
        if(jvalue/=0) &   !20230709
            Fvec(ivalue)=Fvec(ivalue)+Value_observ(jvalue)%value_computation
        if(jvalue==0) & !20230709
            write(7,*)ivalue,Value_observ(ivalue)%value_measure,Value_observ(ivalue)%value_computation
        if(jvalue/=0) & !20230709
            write(7,*)ivalue,Value_observ(ivalue)%value_measure,  &
            (Value_observ(ivalue)%value_computation-Value_observ(jvalue)%value_computation)
    end do
    !   write(7,*)'value_observe%dudx='
    ! do i=1,mvalue
    !write(7,*)i,Value_observ(i)%dudx(:)
    !enddo

    do i=1,n
        do j=1,mvalue
            Fjac(j,i)=Value_observ(j)%dudx(i)
        enddo
    end do
    funx(iter_tr)=dot_product(fvec,fvec)
    do i=1,n
        gk(i)=-2.*dot_product(fvec,Fjac(:,i))
    end do
    !write(7,*)'fvec=',fvec
    !write(7,*)'fjac(:,1)=',fjac(:,1)
    !write(7,*)'fjac(:,2)=',fjac(:,2)
    !write(7,*)'gk=',gk,'funx(iter_tr)=',funx(iter_tr)
20  x1=dot_product(gk,gk)
    x1=sqrt(x1)
    !write(7,*)'x1=',x1
    if(sqrt(x1)<eps) goto 30
    call solve_dx(deltak,N,gk,bk,Sk)
    !write(7,*)'sk=',sk
    if(iter_tr>1)then
        funx(iter_tr)=funx(iter_tr-1)
        gk0=gk
    end if

    x=x0+sk
    do i=1,n
        if(para_back(i)%mode_transform==0)then
            xvalue(i)=x(i)*para_back(i)%factor
        elseif(para_back(i)%mode_transform==1)then
            xvalue(i)=para_back(i)%factor/x(i)
        endif
    end do
    ttime=0.
    Value_observ(:)%value_computation=0.
    call process_analysis
    fvec=0.
    do ivalue=1,m
        if(Value_observ(ivalue)%ic==0)cycle
        Fvec(ivalue)=Value_observ(ivalue)%value_measure-Value_observ(ivalue)%value_computation
        jvalue=Value_observ(ivalue)%jvalue
        if(jvalue/=0) &   !20230709
            Fvec(ivalue)=Fvec(ivalue)-+Value_observ(jvalue)%value_computation
    end do
    bs=Bk.x.sk
    dkstar=funx(iter_tr)-dot_product(fvec,fvec)
    dkpre =-dot_product(gk,sk)-.5*dot_product(sk,bs)

    rk=dkstar/dkpre
    xnew=x0
    !funx(iter_tr)=funx0(iter_tr)
    gk=gk0
    !write(7,*)'dkpre=',dkpre,'dkstar=',dkstar,'rk=',rk
    if(rk>=eta1)then
        xnew=x0+Sk
        do i=1,n
            do J=1,Mvalue
                Fjac(j,i)=Value_observ(j)%dudx(i)
            enddo
        end do
        funx(iter_tr)=dot_product(fvec,fvec)
        if(sqrt(funx(iter_tr))<eps) goto 30
        do i=1,n
            gk(i)=-2.*dot_product(fvec,Fjac(:,i))
        end do
    endif
    !write(7,*)'dkstar=',dkstar,'dkpre=',dkpre,'rk=',rk,'eta1=',eta1,'eta2=',eta2

    !!!校正信赖区间
    if(rk<eta1)then
        deltak=.5*(0+gama1*deltak)
    elseif(rk>=eta1.and.rk<=eta2)then
        deltak=.5*(gama1*deltak+deltak)
    elseif(rk>eta2)then
        if(gama2*deltak<deltab)deltak=.5*(gama2*deltak+deltak)
        if(gama2*deltak>=deltab)deltak=.5*(deltab+deltak)
    endif
    !write(7,*)'deltak=',deltak
    if(rk>=eta1)then   !20230423
        bk0=bk
        yk=gk-gk0
        call update_bk(N,yk,sk,bk0,bk)
    endif
    !!!

    !write(7,*)'x0=',x,'xnew=',xnew
    !!!

    if(iter_tr<mtter)then
        iter_tr=iter_tr+1
        if(iter_tr>=2)then
            !write(7,*)'bk='
            !write(7,*)bk(1,:)
            !write(7,*)bk(2,:)
        endif
        x0=xnew
        goto 20
    end if
    deallocate(gk,gk1,Bk,xnew,bs)
    deallocate(funx,Qs,gk0,bk0,Sk,yk,funx0)

30  continue
    write(7,*)'iter_tr=',iter_tr,'mtter=',mtter
    write(7,*)'abs(df/dx)=',x1,'abs(fvec)=',sqrt(funx(iter_tr))
    write(7,*)'xnew=',xnew

    allocate(fjac1(npoints_pbx,n),errn(n,mobstimes),fjac2(n,n),errx(n,mobstimes))

    allocate(errm0(npoints_pbx,mobstimes),sigma(n))

    errm0=0.
    !write(7,*)'m=',m
    do ivalue=1,m
        if(Value_observ(ivalue)%ic==0)cycle
        iblks_bp=Value_observ(ivalue)%iblks
        iincs_bp=Value_observ(ivalue)%iincs
        istep_bp=Value_observ(ivalue)%istep
        iobse_bp=Value_observ(ivalue)%iobse
        jtotv=Value_observ(ivalue)%ivalue_point
        iobstimes=para_block(iblks_bp)%para_nincs(iincs_bp)%tstep_bp(istep_bp,iobse_bp)
        errm0(jtotv,iobstimes)=fvec(ivalue)
    end do
    !write(7,*)'errm0='
    !      do i=1,mobstimes
    !     write(7,10)errm0(:,i)
    !     end do

10  format(10e15.5)

    do i=1,mobstimes
        fjac1=0.
        do ivalue=1,m
            if(Value_observ(ivalue)%ic==0)cycle
            jtotv=Value_observ(ivalue)%ivalue_point
            iblks_bp=Value_observ(ivalue)%iblks
            iincs_bp=Value_observ(ivalue)%iincs
            istep_bp=Value_observ(ivalue)%istep
            iobse_bp=Value_observ(ivalue)%iobse
            iobstimes=para_block(iblks_bp)%para_nincs(iincs_bp)%tstep_bp(istep_bp,iobse_bp)
            if(i==iobstimes)then
                fjac1(jtotv,:)=fjac(ivalue,:)
            endif
        end do
        !write(7,*)'iobstimes=',i,'fjac1='
        !do jtotv=1,npoints_pbx
        !write(7,10)fjac1(jtotv,:)
        !end do

        do jtotv=1,npoints_pbx
            errn(:,i)=transpose(fjac1).x.errm0(:,i)
        end do
        fjac2=transpose(fjac1).x.fjac1
        allocate(errn0(n,1),errx0(n,1))
        errn0(:,1)=errn(:,i)
        call householderx(fjac2,errn0,errx0)
        errn(:,i)=errn0(:,1)
        errx(:,i)=errx0(:,1)
        deallocate(errn0,errx0)
    end do

    !write(7,*)'errx='
    !    do i=1,mobstimes
    !   write(7,10)errx(:,i)
    !   end do
    !write(7,*)'errn='
    !    do i=1,mobstimes
    !   write(7,10)errn(:,i)
    !   end do
    !

    do j=1,n
        sigma(j)=0.

        do i=1,mobstimes
            sigma(j)=sigma(j)+errx(j,i)**2
        end do
        sigma(j)=sigma(j)/mobstimes
        sigma(j)=sqrt(sigma(j))
    end do


    write(7,*)'sigma=',sigma
    write(7,*)'mobstimes=',mobstimes
    write(7,*)'x=',x

    deallocate(errn,fjac2,fjac1,errx,sigma,errm0)

    end  subroutine trust_region_back_analysis !20220108

    subroutine solve_dx(deltak,N,gk,bk,S)
    integer(ink) N
    real(irk)   a0,b0,c0,yx,lamda1,lamda2,lamda,x1,x2,alfak,gk(:),bk(:,:),bgk(N),  &
        S1(N),Gk1(N,1),S2(N),S21(N,1),skc,skn,deltak,S(:)
    x1=dot_product(gk,gk)
    Bgk=Bk.x.gk
    x2=dot_product(gk,Bgk)
    alfak=x1/x2
    S1=-alfak*gk

    Gk1(:,1)=-Gk
    call householderx(Bk,Gk1,S21)
    S2=S21(:,1)
    skc=dot_product(S1,S1)
    skn=dot_product(S2,S2)
    skc=sqrt(skc)
    skn=sqrt(skn)

    print *,'sqrt(x1)=',sqrt(x1)
    !if(sqrt(x1)<eps) goto 30
    !
    !if(iter_tr==1)then
    !   delta0=x1/10.
    !   write(7,*)'delta0=',delta0,'deltak=',deltak
    !   if(delta0>deltak)deltak=delta0
    !endif
    write(7,*)'s1=',s1
    write(7,*)'s2=',s2
    write(7,*)'deltak=',deltak,'skc=',skc,'skn=',skn
    !write(7,*)'xnew0=',x-deltak*gk/sqrt(x1)
    if(skc>=deltak)then
        S=-deltak*gk/sqrt(x1)
    elseif(skc<deltak.and.skn<=deltak)then
        S=S2   !-S2  20230423
    elseif(skc<deltak.and.skn>deltak)then
        lamda=0.
        a0=dot_product(S2-S1,S2-S1)
        b0=2*dot_product(S2-S1,S1)
        C0=dot_product(S1,S1)
        C0=C0-deltak
        yx=b0**2-4*a0*c0
        lamda1=(-b0+sqrt(yx))/(2*a0)
        lamda2=(-b0-sqrt(yx))/(2*a0)
        write(7,*)'lamda1=',lamda1,'lamda2=',lamda2
        if(lamda1>0.and.lamda2>0.)then
            lamda=lamda1
            if(lamda2>lamda1)lamda=lamda1
        elseif(lamda1>0.)then
            lamda=lamda1
        elseif(lamda2>0.)then
            lamda=lamda2
        endif
        S=S1+lamda*(S2-S1)
    endif

    end subroutine solve_dx

    subroutine update_bk(N,yk,sk,bk0,bk)  !BFGS公式(来源于最优化方法.ppt)
    integer(ink) i,j,N
    real(irk) yk(:),sk(:),bk0(:,:),bk(:,:),BS(N),BS1(N),x1

    !yk=gk-gk0
    !bk0=bk
    Bk=bk0
    bs=Bk0.x.sk

    x1=dot_product(sk,yk)
    do i=1,n
        do j=1,n
            bk(i,j)=bk(i,j)+yk(i)*yk(j)/x1
        end do
    end do

    do i=1,n
        bs1(i)=dot_product(sk,bk0(:,i))
    end do
    x1=dot_product(sk,bs)
    do i=1,n
        do j=1,n
            bk(i,j)=bk(i,j)-bs(i)*bs1(j)/x1
        end do
    end do
    end subroutine update_bk

    subroutine rigid_dis_back_analysis !20211121
    implicit none
    integer(ink) nincs_pb,nstep_pb,mvalue_point,imt,igdis_bk,i,ipoin,istep,npara,nmbpoint
    real   (irk),allocatable:: observ(:),observt(:,:),interp(:,:),interpt(:,:),rgdis(:,:)
    integer(ink),allocatable::measurep(:),backrdisp(:),changepx(:)
    real   (irk) ug,ue
    allocate(measurep(Npoints_pbx))
    do iblks=1,Nblks_pb
        nincs_pb=para_block(iblks)%nincs_pb
        do iincs=1,nincs_pb
            nstep_pb=para_block(iblks)%para_nincs(iincs)%nstep_pb
            mvalue_point=para_block(iblks)%para_nincs(iincs)%mvalue_point

            measurep=0
            do imt=1,mvalue_point
                ipoin=para_block(iblks)%para_nincs(iincs)%list_point_pb(imt)
                measurep(ipoin)=1
            end do

            allocate(backrdisp(Npoints_pbx),changepx(Npoints_pbx))

            do igdis_bk=1,ngdis_bk
                write(7,*)'iblks=',iblks,'iincs=',iincs,'igdis_bk=',igdis_bk
                write(7,*)'      istep                 刚体位移'

                backrdisp=0
                backrdisp(rigid_bk(igdis_bk)%node_bk)=1

                nmbpoint=0
                changepx=0

                do ipoin=1,Npoints_pbx
                    if(backrdisp(ipoin)==1.and.measurep(ipoin)==1)then
                        nmbpoint=nmbpoint+1
                        changepx(ipoin)=nmbpoint
                    endif
                end do

                npara=3*(ndimn-1)
                allocate(observ(nmbpoint),observt(npara,1),interp(nmbpoint,npara),interpt(npara,npara),rgdis(npara,1))

                do istep=1,nstep_pb !istep
                    nmbpoint=0
                    do i=1,mvalue  !i
                        if(iblks/=Value_observ(i)%iblks)cycle
                        if(iincs/=Value_observ(i)%iincs)cycle
                        if(istep/=Value_observ(i)%istep)cycle
                        ipoin=Value_observ(i)%ivalue_point
                        if(backrdisp(ipoin)/=1.or.measurep(ipoin)/=1)cycle
                        nmbpoint=nmbpoint+1
                        observ(nmbpoint)=Value_observ(i)%value_measure
                        idofn=Value_observ(i)%idofn   !20230523
                        interp(nmbpoint,1:npara)=rigid_bk(igdis_bk)%npdisp(idofn,changepx(ipoin),:)
                    end do !i
                    interpt=transpose(interp).x.interp
                    observt(:,1)=transpose(interp).x.observ
                    do idofn=1,npara
                        if(rigid_bk(igdis_bk)%fixed_dis(idofn)==0) then
                            if(interpt(idofn,idofn)<1.e-1)then
                                interpt(idofn,idofn)=1.e10
                            else
                                interpt(idofn,idofn)=interpt(idofn,idofn)*1.e10
                            endif
                        endif
                    end do
                    call householder(interpt,observt,rgdis)


                    write(7,20)istep,rgdis(:,1)
                    write(7,*)'    观测点号    刚体位移      弹性位移'
                    nmbpoint=0
                    do i=1,mvalue
                        if(iblks/=Value_observ(i)%iblks)cycle
                        if(iincs/=Value_observ(i)%iincs)cycle
                        if(istep/=Value_observ(i)%istep)cycle
                        ipoin=Value_observ(i)%ivalue_point
                        if(backrdisp(ipoin)/=1.or.measurep(ipoin)/=1)cycle
                        nmbpoint=nmbpoint+1
                        ug=dot_product(interp(nmbpoint,:),rgdis(:,1))
                        ue=observ(nmbpoint)-ug
                        write(7,20)i,ug,ue
                    end do


                    !        do jpoin=1,npara
                    !ipoin=nodvar_bk(igdis_bk)%node_bk(jpoin)
                    !itotv=nodfn(lmdofn(kdofn),ipoin)
                    !result_zero(itotv)=nodvar(jpoin,1)
                    !end do
                    if(outplot(1:3)=='GID')   call OUT_GID_WRITE

                end do !istep
                deallocate(observ,observt,interp,interpt,rgdis,changepx)
            end do  !igdis_bk
            deallocate(backrdisp)
        end do !iincs
    end do !iblks
    deallocate(measurep)

20  format(i10,6e15.5)

    end subroutine rigid_dis_back_analysis !20211121

    subroutine nodal_value_back_analysis !20211201
    implicit none
    integer(ink) nincs_pb,nstep_pb,mvalue_point,imt,igdis_bk,i,ipoin,istep,   &
        npara,nmbpoint,jdofn,itotv,kdofn,kdofn1
    real   (irk),allocatable:: observ(:),observt(:,:),interp(:,:),interpt(:,:),nodvar(:,:)
    integer(ink),allocatable::measurep(:),backrdisp(:),listnode(:)

    allocate(measurep(Npoints_pbx))
    do iblks=1,Nblks_pb
        nincs_pb=para_block(iblks)%nincs_pb
        do iincs=1,nincs_pb
            nstep_pb=para_block(iblks)%para_nincs(iincs)%nstep_pb
            mvalue_point=para_block(iblks)%para_nincs(iincs)%mvalue_point

            measurep=0
            do imt=1,mvalue_point
                ipoin=para_block(iblks)%para_nincs(iincs)%list_point_pb(imt)
                measurep(ipoin)=1
            end do

            allocate(backrdisp(Npoints_pbx),listnode(npoin))

            do igdis_bk=1,ngval_bk
                write(7,*)'iblks=',iblks,'iincs=',iincs,'igvar_bk=',igdis_bk
                write(7,*)'      istep','      node_value'

                backrdisp=0
                listnode=0


                listnode(nodvar_bk(igdis_bk)%node_bk(:))=1

                do i=1,Npoints_pbx
                    nintf=para_points(i)%nintf
                    if(nintf==0)cycle
                    if(all(listnode(para_points(i)%listf(:))==1))backrdisp(i)=1
                end do

                nmbpoint=0

                do ipoin=1,Npoints_pbx
                    if(backrdisp(ipoin)==1.and.measurep(ipoin)==1)nmbpoint=nmbpoint+1
                end do

                npara=nodvar_bk(igdis_bk)%npoin_bk
                allocate(observ(nmbpoint),observt(npara,1),interp(nmbpoint,npara),interpt(npara,npara),nodvar(npara,1))
                interp=0.;interpt=0.
                do istep=1,nstep_pb !istep

                    do kdofn1=1,nodvar_bk(igdis_bk)%ndofn_bk
                        kdofn=nodvar_bk(igdis_bk)%dof_bk(kdofn1)

                        nmbpoint=0
                        do i=1,mvalue  !i
                            if(iblks/=Value_observ(i)%iblks)cycle
                            if(iincs/=Value_observ(i)%iincs)cycle
                            if(istep/=Value_observ(i)%istep)cycle
                            if(kdofn/=Value_observ(i)%idofn)cycle
                            ipoin=Value_observ(i)%ivalue_point
                            if(backrdisp(ipoin)/=1.or.measurep(ipoin)/=1)cycle
                            nmbpoint=nmbpoint+1
                            observ(nmbpoint)=Value_observ(i)%value_measure
                            do  idofn=1,para_points(ipoin)%nintf
                                jpoin=para_points(ipoin)%listf(idofn)
                                jdofn=nodvar_bk(igdis_bk)%nodet_bk(jpoin)
                                interp(nmbpoint,jdofn)=para_points(ipoin)%rintf(idofn)
                            end do
                        end do !i


                        interpt=transpose(interp).x.interp
                        observt(:,1)=transpose(interp).x.observ
                        call householder(interpt,observt,nodvar)
                        write(7,20)istep,nodvar(:,1)
                        do jpoin=1,npara
                            ipoin=nodvar_bk(igdis_bk)%node_bk(jpoin)
                            itotv=nodfn(lmdofn(kdofn),ipoin)
                            result_zero(itotv)=nodvar(jpoin,1)
                        end do


                        if(outplot(1:3)=='GID')   call OUT_GID_WRITE
                    end do !kdofn1
                end do !istep
                deallocate(observ,observt,interp,interpt,nodvar)
            end do  !igdis_bk
            deallocate(backrdisp)
        end do !iincs
    end do !iblks
    deallocate(measurep,listnode)

10  format(20e15.3)
20  format(i10,20e15.3)

    end subroutine nodal_value_back_analysis !20211201

    subroutine parameter_back_analysis_verify !20200812
    implicit none

    INTEGER   k0,nincs_pb,i,j,istep,inode,jnode,idofn,mvalue_point, &
        iblks,iincs,nintf,mode_transform,bblks,i0,ibpstep,mback_point, &
        nstep_pb,ix
    real     factor,dtime_pb
    integer,allocatable::list_point_pb(:)

    do istoch=1,nstoch
        xvalue=para_stoch(:,istoch)
        do i=1,npara
            factor=para_back(i)%factor
            mode_transform=para_back(i)%mode_transform
            if(mode_transform==0)then
                xvalue(i)=factor*xvalue(i)
            elseif(mode_transform==1)then
                xvalue(i)=factor/xvalue(i)
            endif
        end do
        tbstep=0
        ttime=0.
        print *,'istoch=',istoch
        call process_analysis
    end do



    !write(observ_unit,*)'information for given points：Npoints_pb'
    !write(observ_unit,10) nback_point
    !write(observ_unit,*)'1:Npoints_pb/i0,ndofn,imdofn,nintf'
    ! k0=0 ; nintf=1
    !              do i=1,nback_point
    !        k0=k0+1
    !        inode=freedom_for_back(1,i)
    !        idofn=freedom_for_back(2,i)
    !           write(observ_unit,10)k0,1,idofn,1
    !           write(observ_unit,10)inode
    !           write(observ_unit,12)1.
    !              end do
    write(observ_unit,*)'Nblks_pb/1:nblks_pb->/text/nincs_pb/1:nincs_pb->dtime_pb,nstep_p,bobserv_pb'
    write(observ_unit,10)runblks
    rewind(mainunit)
    do iblks=1,runblks
        read(mainunit,*)text
        write(observ_unit,*)text
        read(mainunit,*)nincs_pb
        write(observ_unit,10)nincs_pb
        do iincs=1,nincs_pb
            read(mainunit,*)i0,dtime_pb,i0,i0,nstep_pb
            write(observ_unit,13)dtime_pb,nstep_pb,1,nstoch
            read(mainunit,*)text
        end do
    end do

    write(observ_unit,*)'observed values(1:nback_point)|1:nstoch:(1:tbstep)'
    ix=0
    do i=1,nback_point
        do j=1,nstoch
            ix=ix+1
            write(observ_unit,10)ix,i,j,0
            write(observ_unit,12)Value_vc(i,:,j)
        end do
    end do

10  format(20i10)
12  format(20e15.5)
13  format(e15.5,3i10)

    end  subroutine parameter_back_analysis_verify

    subroutine  parameter_back_analysis_read   !20190810

    implicit none
    integer             I, J,k,nincs_pb,nstep_pb,ivalue_point,i0,inode,idofn,jnode,nintf
    INTEGER             nvalue,mvalue_point,k0,observ_pb,iobse,ndofn,ivalue,i1,j1,j0,ibstep
    INTEGER,allocatable::listdofn(:)
    real(irk)           dtime_pb,ttime_pb
    real(irk),allocatable::obs_value(:)


    read(back_ctl_unit,*)text
    read(back_ctl_unit,*)Npara  !20230523


    allocate(Xvalue(Npara),Para_back(Npara))  !20230523
    read(back_ctl_unit,*)text
    do i=1,Npara
        read(back_ctl_unit,*)j,para_back(i)%imat,para_back(i)%name,para_back(i)%factor,  &
            para_back(i)%factor_inc,para_back(i)%mode_transform
    enddo
    read(back_ctl_unit,*)text
    read(back_ctl_unit,*)eps,iter1,iter2,rs,jac_eps
    read(back_ctl_unit,*)text  !待反演参数初始值
    read(back_ctl_unit,*)Xvalue

    write(7,*)'xvalue=',xvalue
    write(7,*)'eps,iter1,iter2,rs,jac_eps=',eps,iter1,iter2,rs,jac_eps

    end subroutine parameter_back_analysis_read  !20190810

    subroutine  trust_region_back_analysis_read  !20220108

    implicit none
    integer             I, J,k,nincs_pb,nstep_pb,ivalue_point,i0,inode,idofn,jnode,nintf
    INTEGER             nvalue,mvalue_point,k0,observ_pb,iobse,ndofn,ivalue,j0,i1,j1,ibstep
    INTEGER,allocatable::listdofn(:)
    real(irk)           dtime_pb,ttime_pb
    real(irk),allocatable::obs_value(:)


    allocate(trustp(1))
    read(back_ctl_unit,*)text
    read(back_ctl_unit,*)Npara  !20230523


    allocate(Xvalue(Npara),Para_back(Npara))  !20230523
    read(back_ctl_unit,*)text

    do i=1,Npara
        read(back_ctl_unit,*)j,para_back(i)%imat,para_back(i)%name,para_back(i)%factor,  &
            para_back(i)%factor_inc,para_back(i)%mode_transform
    enddo

    read(back_ctl_unit,*)text
    print *,text
    read(back_ctl_unit,*)trustp(1)%eta1,trustp(1)%eta2,   &
        trustp(1)%gama1,trustp(1)%gama2,trustp(1)%eps,  &
        trustp(1)%eta01,trustp(1)%eta02,trustp(1)%delta0, &
        trustp(1)%deltab,trustp(1)%mtter
    read(back_ctl_unit,*)text  !待反演参数初始值
    read(back_ctl_unit,*)Xvalue


    end subroutine trust_region_back_analysis_read  !20220108

    subroutine  parameter_back_analysis_verify_read   !20200813

    implicit none
    integer             I, J,k0,k,tbstep,nincs_pb,nstep_pb,bblks
    !INTEGER             nvalue,mvalue_point,k0
    !DOUBLE PRECISION    value
    real    timebprecord,dtime_pb
    real,allocatable:: ruo(:,:),sigma_t(:)

    open(back_ctl_unit,file=probn(1:len1)//'.btl')
    open(observ_unit,file=probn(1:len1)//'.obsc')

    read(back_ctl_unit,*)text  !输入与Bparamete为负时的相关内容（给定随机变量参数，进行正分析计算，输出相关结果用于反演分析方法验证）

    read(back_ctl_unit,*)Npara,nback_point,nstoch
    print *,'Npara,nback_point,nstoch=',Npara,nback_point,nstoch
    tbstep=0
    nblks_pb=runblks
    allocate(para_block(Nblks_pb))
    rewind(mainunit)
    do i=1,nblks_pb
        read(mainunit,*)text
        read(mainunit,*)nincs_pb
        para_block(i)%nincs_pb=nincs_pb
        allocate(para_block(i)%para_nincs(nincs_pb))
        timebprecord=0.
        do j=1,nincs_pb
            read(mainunit,*)k0,dtime_pb,k0,k0,nstep_pb
            read(mainunit,*)text
            para_block(i)%para_nincs(j)%nstep_pb=nstep_pb
            tbstep=tbstep+nstep_pb
            allocate(para_block(i)%para_nincs(j)%time_bp(nstep_pb))
            do k=1,nstep_pb
                timebprecord=timebprecord+dtime_pb
                para_block(i)%para_nincs(j)%time_bp(k)=timebprecord
            end do
        end do
    end do


    allocate(Para_back(Npara),freedom_for_back(4,nback_point))
    allocate(Xvalue(Npara),mean_value(Npara),sigma_value(Npara),para_stoch(npara,nstoch))

    read(back_ctl_unit,*)text
    print *,text
    do i=1,Npara
        read(back_ctl_unit,*)j,para_back(i)%imat,para_back(i)%name,para_back(i)%factor,para_back(i)%mode_transform
    enddo

    read(back_ctl_unit,*)text
    read(back_ctl_unit,*)mean_value  !随机变量均值
    read(back_ctl_unit,*)sigma_value !随机变量离散系数
    do i=1,Npara
        sigma_value(i)=sigma_value(i)*mean_value(i)
    end do
    if(nstoch==1)then !20210805
        para_stoch(:,1)=mean_value
    else !20210805
        do i=1,Npara
            call vsl_gauss_gen_single(nstoch,mean_value(i),sigma_value(i),para_stoch(i,:))
        end do
    endif !20210805

    write(7,*)'蒙特卡洛分析随机抽样参数生成结果'
    do i=1,nstoch
        write(7,10)i,para_stoch(1:Npara,i)
    end do

    allocate(ruo(npara,npara),sigma_t(npara))
    ruo=0.
    do i=1,npara
        xvalue(i)=sum(para_stoch(i,:))
        xvalue(i)=xvalue(i)/nstoch
    end do

    sigma_t=0.
    do i=1,npara
        do j=1,nstoch
            sigma_t(i)=sigma_t(i)+(para_stoch(i,j)-xvalue(i))**2
        end do
    end do

    !write(7,*)'sigma_t=',sigma_t
    if(nstoch>1)then  !20210805
        do i=1,npara
            do j=1,npara
                do k=1,nstoch
                    ruo(i,j)=ruo(i,j)+(para_stoch(i,k)-xvalue(i))*(para_stoch(j,k)-xvalue(j))
                    !if(i==2.and.j==3)write(7,*)'ruo(i,j)=',ruo(i,j)
                end do
                ruo(i,j)= ruo(i,j)/sqrt(sigma_t(i)*sigma_t(j))
            end do
        end do

        do i=1,npara
            sigma_t(i)=sigma_t(i)/(nstoch-1)
            sigma_t(i)=sqrt(sigma_t(i))
        end do
    endif !20210805


    write(7,*)'模拟数据均值与与给定值比较'
    write(7,11)mean_value  !随机变量均值
    write(7,11)xvalue
    write(7,*)'模拟数据方差与与给定值比较'
    write(7,11)sigma_value  !随机变量均值
    write(7,11)sigma_t
    write(7,*)'相关系数'
    do i=1,npara
        write(7,11)ruo(i,:)
    end do


    deallocate(ruo,sigma_t)


    read(back_ctl_unit,*)text
    do i=1,nback_point
        read(back_ctl_unit,*)k0,inode,idofn,jnode,bblks
        freedom_for_back(1,i)=inode
        freedom_for_back(2,i)=idofn
        freedom_for_back(3,i)=jnode
        freedom_for_back(4,i)=bblks
    end do

    write(7,*)'nback_point=',nback_point,'tbstep,nstoch=',tbstep,nstoch
    allocate(value_vc(nback_point,tbstep,nstoch))
    value_vc=99999.
    !write(7,*)'value_vc(nback_point,tbstep,nstoch)=',value_vc(nback_point,tbstep,nstoch)
10  format(i10,20e15.5)
11  format(20e15.5)

    end subroutine parameter_back_analysis_verify_read   !20200813

    subroutine  observe_back_analysis_read  !20230523

    implicit none
    integer             I, J,k,nincs_pb,nstep_pb,ivalue_point,i0,j0,inode,idofn,jnode,nintf
    INTEGER             nvalue,mvalue_point,k0,observ_pb,iobse,ndofn,ivalue,i1,j1,ibstep,  &
        ix,jvalue,begin_day_obs,dstep_pb,dtime_pb, &
        begin_day_pb,end_day_pb,ttime_pb   !20230709
    INTEGER,allocatable::listdofn(:)
    DOUBLE PRECISION    value

    real(irk),allocatable::obs_value(:)
    ! ttime_pb,dtime_pb的单位为天，所以实际工程反分析时将这两个变量按整数处理。

    !open(observ_unit,file=probn(1:len1)//'.obs')
    !read(back_ctl_unit,*)text  !输入与Bparameter/=0时的相关内容（不为零时，执行参数优化反演）
    tbstep=0
    ttime_pb=0.
    mobstimes=0
    Mvalue=0
    read(observ_unit,*)text
    read(observ_unit,*)Nblks_pb  !20230523
    allocate(para_block(Nblks_pb))  !20230523

    do iblks=1,Nblks_pb
        read(observ_unit,*)text
        read(observ_unit,*) nincs_pb
        para_block(iblks)%nincs_pb=nincs_pb
        allocate(para_block(iblks)%para_nincs(nincs_pb))
        do i=1,nincs_pb
            read(observ_unit,*) dtime_pb,nstep_pb,observ_pb,begin_day_pb,end_day_pb  !20230924
            dstep_pb=1
            tbstep=tbstep+nstep_pb !20230719
            para_block(iblks)%para_nincs(i)%nstep_pb=nstep_pb
            para_block(iblks)%para_nincs(i)%dtime_pb=dtime_pb
            para_block(iblks)%para_nincs(i)%dstep_pb=dstep_pb
            para_block(iblks)%para_nincs(i)%begin_day_pb=begin_day_pb
            para_block(iblks)%para_nincs(i)%end_day_pb=end_day_pb
            para_block(iblks)%para_nincs(i)%mvalue_point=Npoints_pbx
            mvalue_point=Npoints_pbx
            allocate(para_block(iblks)%para_nincs(i)%time_bp(nstep_pb),  &
                para_block(iblks)%para_nincs(i)%tstep_bp(nstep_pb,observ_pb),  &
                para_block(iblks)%para_nincs(i)%list_point_pb(mvalue_point))

            do j=1,Npoints_pbx
                para_block(iblks)%para_nincs(i)%list_point_pb(j)=j
            end do

            do j=1,nstep_pb,dstep_pb !20230719
                ttime_pb=begin_day_pb+dtime_pb*j*dstep_pb !20230924
                para_block(iblks)%para_nincs(i)%time_bp(j)=ttime_pb
                do i0=1,observ_pb
                    mobstimes=mobstimes+1
                    para_block(iblks)%para_nincs(i)%tstep_bp(j,i0)=mobstimes
                end do
            end do  !j
            do j=1,Npoints_pbx
                Mvalue=mvalue+((nstep_pb-1)/dstep_pb+1)*observ_pb   !20230719
            end do
        end do !i
    end do  !iblks

    allocate(Value_observ(Mvalue))

    read(observ_unit,*)text
    print *,'text=',text
    allocate(obs_value(tbstep))
    ivalue=0

    do j=1,Npoints_pbx

        !begin_day_obs=para_points(j)%begin_day_obs

        do j0=1,observ_pb
            read(observ_unit,*)ix,i1,j1,begin_day_obs    !对应测点点号，方向号
            !print *,'j=','j0=',j0,'i1=',i1,'j1=',j1
            read(observ_unit,*)obs_value

            ibstep=0   !20230902
            do iblks=1,Nblks_pb
                nincs_pb=para_block(iblks)%nincs_pb
                do iincs=1,nincs_pb
                    nstep_pb=para_block(iblks)%para_nincs(iincs)%nstep_pb
                    jvalue=0
                    do istep=1,nstep_pb,dstep_pb  !20230719
                        ibstep=ibstep+1  !20230902
                        ivalue=ivalue+1
                        !if(begin_day_obs==para_block(iblks)%para_nincs(iincs)%time_bp(istep))   &
                        if(begin_day_obs==ibstep) &
                            jvalue=ivalue   !20230924
                        Value_observ(ivalue)%iblks=iblks
                        Value_observ(ivalue)%iincs=iincs
                        Value_observ(ivalue)%istep=istep
                        Value_observ(ivalue)%iobse=j0
                        Value_observ(ivalue)%idofn=para_points(j)%listdofn(1)
                        Value_observ(ivalue)%ivalue_point=j
                        if(Value_observ(ivalue)%idofn<=ndimn)  & !20230709
                            Value_observ(ivalue)%jvalue=jvalue  !20230709

                        Value_observ(ivalue)%ic=1
                        if(abs(obs_value(ibstep)-99999.)<.01)Value_observ(ivalue)%ic=0   !20230902

                        if(Value_observ(ivalue)%idofn<=ndimn.and.(ivalue==jvalue))  &
                            Value_observ(ivalue)%ic=0   !20230709

                        Value_observ(ivalue)%value_measure=obs_value(ibstep)  !20230902


                        if(balgor/=0) then !20220108
                            allocate(Value_observ(ivalue)%dudx(npara))
                            Value_observ(ivalue)%dudx=0.
                        endif
                    end do  !istep
                end do  !iincs
            end do !iblks
        end do !j0
    end do !j

    deallocate(obs_value)
    if(Bparameter==1.or.Bparameter==2)then
        allocate(Fvec(Mvalue),Fvec1(Mvalue),Fvec2(Mvalue),FJAC (Mvalue,Npara))
        fvec=0.
        fjac=0.
    endif
    !

    print *,'mvalue=',mvalue


    end subroutine observe_back_analysis_read  !20230523

    SUBROUTINE back_analysis  !20150925

    logical logx
    character(80)text
    integer(ink) itotv,ielem,irst,trstep0,ipoin,idofn,ij,idofix,ldofix,idelgroup,i0,ipairs,k
    real   (irk) xtime,time_begin,detal,ttime0,coef,xij,qi,qerr,qabs
    real   (irk) dispoint0,dispoint1,tQerr,tQabs,dispoint1g,dispoint0g  !20210726
    real   (irk) dis2e,diser,kmodu

    real   (irk),allocatable::rvectorm(:),value(:)
    integer(ink) iintf,nintf,iieq,igapbf,mdism   !!int2000

    integer(ink) i,iincs_i,iblks_i,istep_i,inode,jnode,ivalue,ivalue_point,bblks  !20230523
    integer(ink),pointer::listf(:)  !20200819
    real   (irk),pointer::rintf(:)  !20200819


    integer(ink) igapb,npgblock,jpoin,igaps,ipair,idimn,itotvbt,jdimn, &  !!ctt2005
        jtotv,kpoin,lpoin,jtotvbt,npairs,cwater,jgaps,jpair,kdimn,jpoin0,itotv0   !!ctt2005
    real   (irk),allocatable::rot(:,:),tofor0(:),uireact(:,:)  !!ctt2005
    real   (irk),allocatable::unitl(:),unitg(:),cmatrixl(:,:) !!ctt2005

    if (meshc==1.or.rmesh/=0)rewind(mainunit)
    if(Bparameter/=0.and.iblks==1)rewind(mainunit)  !20190810

    !read(observ_unit,*)text


    read(mainunit,*)text
    read(mainunit,*)nincs

    print *,' in back_analysis'
    tQerr=0.;tQabs=0.

    if(ngaps/=0.or.nrcsteel/=0)allocate(tofor0(ntotv)) !!ctt2005

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


        ttime0=ttime
        trstep0=trstep
        !if(nbackf/=0) & !20210726,20230430
        !read(mainunit,*)text

        do istep=inc_step,nstep,inc_step

            if(iblks>=stab_matde)call stab_initialize

            write(chkunit,*)'Increment step=',istep
            if(nbackf/=0)trstep=trstep0+istep  !20220101
            if(outintr>0.and.iblks>=outintr)trstep=trstep0+istep  !20200226
            xtime=ditime*istep
            ttime=ttime0+ditime*istep !! only for output

            result_zero=0.0  !201605
            do igapb=1,ngapb
                if(gapb(igapb)%nrdof/=0)gapb(igapb)%rdisp_zero=0.
            end do


            do igaps=1,ngaps   !20190810
                npairs=gaps(igaps)%npairs
                do ipairs=1,npairs
                    if(gaps(igaps)%pair_process(ipairs)==0)cycle  !20200331
                    gaps(igaps)%dxyz0(:,ipairs)=0.
                    gaps(igaps)%dxyz(:,ipairs)=0.
                    gaps(igaps)%ctforce0(:,ipairs)=0.
                    gaps(igaps)%ctforce(:,ipairs)=0.
                end do
            end do   !20190810
            call gpvar_initial  !201605

            if(nbackf/=0)then  !20210804
                do idofn=1,nbackf  !20210804
                    mdism=backf(idofn)%mdism
                    backf(idofn)%ic=1
                    do i0=1,mdism
                        if(abs(backf(idofn)%dism(i0,trstep)-99999.)<.01) &  !20220101
                            backf(idofn)%ic(i0)=0
                    end do
                enddo !20210804
            endif  !20210804


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
                if(type_load=='DISCONTROL')preact0=prescrib(1)%rdofix
                if((ngaps/=0.or.nrcsteel/=0).and.mdiv==1)tofor0=tofor !!ctt2005
                if((ngaps/=0.or.nrcsteel/=0).and.mdiv/=1)tofor0=toform !!ctt2005
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


1                   if(type_load=='LOAD2'.or.(kstat==2.and.iiter.le.2).or.(type_load/='LOAD2'.and.kstat/=2.and.iiter==1).or.  &
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

                        call eload_initialize

                        call residu_f

                        if(rmesh>0.and.nelem1>0)call residu_f1
                        if(rmesh>1.and.nelem2>0)call residu_f2
                        call eload_field
                        if(ground_inf/=0)call semi_inf_load
                        call force_internal
                    endif !for iiter==1 and istep==inc_step .and.idiv==1  temperature


                    if(ngaps/=0.and.iblks>=iblks_bt.and.iiter==1.and.mdiv==1)call ctfor_to_tofor(tofor0,tofor)  !!ctt2005
                    if(ngaps/=0.and.iblks>=iblks_bt.and.iiter==1.and.mdiv/=1)call ctfor_to_tofor(tofor0,toform)  !!ctt2005
                    if(nrcsteel/=0.and.iiter==1.and.mdiv==1)call csfor_to_tofor(tofor0,tofor)  !!20210328
                    if(nrcsteel/=0.and.iiter==1.and.mdiv/=1)call csfor_to_tofor(tofor0,toform)  !!20210328

                    if (mdiv/=1) then
                        if(idiv==1.and.iiter==1.and.allocated(torel))toform=toform+torel  !!20210328
                    else
                        if(iiter==1.and.allocated(torel))tofor=tofor+torel !!20210328
                    endif


                    if(neuman==1.and.(istep/=1.or.iiter/=1).and.kresl/=0)  goto 2  !ctt2005 , change position!
                    if(type_nl==8.and.(iiter>1.or.(kstat==2.and.iiter>2))) goto 2  !MNR
                    if ((type_solver=='PROFILE'.or.type_solver=='PARDISO').and.kresl/=0)then
                        operation='FACTORIZE'
                        call solve
                    end if
2                   continue

                    logx=ngaps/=0.and.(iiter==1.and.istep==inc_step.and.iincs==(lincs+1)).and.iblks==iblks_bt
                    if (logx)then !ctt2005
                        if (restart_ctt==0)then !restart_ctt
                            kdimn=ndimn
                            if(block_stab==1)kdimn=3*(ndimn-1) !2015/11/17
                            allocate(rot(kdimn,kdimn))
                            rot=0.
                            do igapbf=1,nbackf
                                igapb=backf(igapbf)%groupb
                                npgblock=gapb(igapb)%npgblock
                                if(backf(igapbf)%mdism>npgblock*kdimn) allocate(uireact(backf(igapbf)%mdism,npgblock*kdimn))
                                gapb(igapb)%cmatrix=0.
                                do ipoin=1,npgblock
                                    igaps=gapb(igapb)%nodegblock_igaps(ipoin)
                                    ipair=gapb(igapb)%nodegblock_ipairs(ipoin)
                                    ij=gapb(igapb)%nodegblock_onetwo(ipoin)
                                    coef=1.
                                    if(ij==2)coef=-1.
                                    rot(1:ndimn,1:ndimn)=gaps(igaps)%rot(:,:,ipair)  !2015/11/17
                                    if(kdimn>ndimn)then  !2015/11/17
                                        if(ndimn==2)rot(3,3)=1.
                                        if(ndimn==3)rot(4:6,4:6)= rot(1:ndimn,1:ndimn)
                                    endif
                                    allocate(unitl(kdimn),unitg(kdimn))
                                    do idimn=1,kdimn
                                        itotvbt=(ipoin-1)*kdimn+idimn
                                        unitl=0.
                                        unitl(idimn)=1.*coef
                                        unitg=transpose(rot).x.unitl  !2015/11/17
                                        rvector=0.
                                        call  unit_force_trans(igapbf,kdimn,ij,unitg,igaps,ipair,rvector)

                                        operation='SOLVE'
                                        call solve

                                        if(backf(igapbf)%mdism==npgblock*kdimn)then
                                            do kpoin=1,backf(igapbf)%mdism
                                                !jdimn=backf(igapbf)%listdim(kpoin)
                                                dispoint1=0.    !20210726
                                                nintf=backf(igapbf)%relat(kpoin)%nintf

                                                do iintf=1,nintf
                                                    jtotv=backf(igapbf)%relat(kpoin)%listf(iintf)
                                                    dispoint1=dispoint1+result(jtotv)*backf(igapbf)%relat(kpoin)%rintf(iintf)
                                                end do
                                                jtotvbt=kpoin
                                                gapb(igapb)%cmatrix(jtotvbt,itotvbt)=dispoint1
                                            end do  !kpoin
                                        else if(backf(igapbf)%mdism>npgblock*kdimn)then
                                            do kpoin=1,backf(igapbf)%mdism
                                                !jdimn=backf(igapbf)%listdim(kpoin)
                                                !itotv=nodfn(jdimn,jpoin)
                                                dispoint1=0.    !20210726
                                                nintf=backf(igapbf)%relat(kpoin)%nintf
                                                do iintf=1,nintf
                                                    jtotv=backf(igapbf)%relat(kpoin)%listf(iintf)
                                                    dispoint1=dispoint1+result(jtotv)*backf(igapbf)%relat(kpoin)%rintf(iintf)
                                                end do

                                                jtotvbt=kpoin
                                                uireact(jtotvbt,itotvbt)=dispoint1
                                            end do  !kpoin

                                        endif
                                    end do  !idimn
                                    deallocate(unitl,unitg)
                                end do  !ipoin

                                if(backf(igapbf)%mdism>npgblock*kdimn)then
                                    gapb(igapb)%uireact=uireact
                                    gapb(igapb)%cmatrix(1:npgblock*kdimn,1:npgblock*kdimn)=transpose(uireact).x.gapb(igapb)%uireact
                                    deallocate(uireact)
                                endif


                            end do  !igapb
                            call forAdirect_back_analysis !fzx !形成A矩阵


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
                            call forAdirect_back_analysis !fzx !形成A矩阵
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
                    rvector=0.0
                    !write(7,*)'iiter=',iiter,'tofor,stfor,rvector='
                    if (type_solver/='JPCG') then
                        do itotv=1,ntotv
                            if (totveq(itotv)/=0)then
                                if (mdiv/=1)then
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                        toform(itotv)-stfor(itotv)

                                else
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                        tofor(itotv)-stfor(itotv)

                                    !write(7,*)itotv,tofor(itotv),stfor(itotv), rvector(totveq(itotv))
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


                    if (type_nl==8)then
                        if(kstat/=2)call bfgsr(iiter)
                        if(kstat==2)call bfgsr(iiter-1)
                    else
                        operation='SOLVE'
                        call solve
                    endif
                    result_zero_e=result_zero_e+result

                    if(ngaps/=0.and.iblks>=abs(iblks_bt))call solve_ctt_back_analysis  !!ctt2005

                    !write(7,*)'result='
                    !do itotv=1,ntotv
                    !write(7,*)itotv,result(itotv)
                    !end do


                    call TIME(char_time)
                    print *, 'time: ', char_time
                    write(chkunit,*)'time: ', char_time


                    call varupdate
                    call eload_initialize
                    call residu_f
                    if(rmesh>0.and.nelem1>0)call residu_f1
                    if(rmesh>1.and.nelem2>0)call residu_f2
                    if(type_load/='LOAD2'.or.(type_load=='LOAD2'.and.idiv==2))then   !806
                        call eload_field
                        if(ground_inf/=0)call semi_inf_load

                        call reaction_prescribed


                        call conver_load
                        if(nchek==0) call conver_nodal_value

                        if(nchek==0)exit !tcl
                    endif !ep2010
10                  continue
                    print *,'miter=',miter,'iiter=',iiter
                end do   !! loop for iiter

                if(type_load/='LOAD2') &
                    call gpvarupdate

                if(rmesh>0.and.nelem1>0)call gpvarupdate1
                if(rmesh>1.and.nelem2>0)call gpvarupdate2
            end do    !! for idiv
            if(type_load=='LOAD2') &
                call gpvarupdate
            if(modf_dis_blocks(iblks)==1)call construction_dis_modify


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

            do igapbf=1,nbackf
                write(7,*)'nodal number,direction,ratio,rigid_dis,elastic_dis'
                dis2e=0.;diser=0.
                qi=0.
                do kpoin=1,backf(igapbf)%mdism
                    if(backf(igapbf)%ic(kpoin)==0)cycle  !20230523
                    jdimn=backf(igapbf)%listdim(kpoin)
                    dispoint1=0.;dispoint1g=0.
                    nintf=backf(igapbf)%relat(kpoin)%nintf
                    do iintf=1,nintf
                        jtotv=backf(igapbf)%relat(kpoin)%listf(iintf)
                        dispoint1=dispoint1+result_zero_e(jtotv)*backf(igapbf)%relat(kpoin)%rintf(iintf)
                        dispoint1g=dispoint1g+(result_zero(jtotv)-result_zero_e(jtotv))*backf(igapbf)%relat(kpoin)%rintf(iintf)

                    end do

                    qi=abs(dispoint1g/dispoint1)
                    !dis2e=dis2e+(dispoint1-dispoint1g)**2
                    !diser=diser+(backf(igapbf)%dism(kpoin,trstep)-dispoint1g)*(dispoint1-dispoint1g)
                    write(7,20)kpoin,jdimn,qi,dispoint1g,dispoint1
                end do
            end do

            !    kmodu=diser/dis2e
            !write(7,*)'kmodu=',kmodu

            do igapbf=1,nbackf
                write(7,*)'relative displacement error/nodal number,direction,relative error,u_observation,u_computation'
                qi=0.; Qerr=0.;Qabs=0.
                do kpoin=1,backf(igapbf)%mdism
                    if(backf(igapbf)%ic(kpoin)==0)cycle  !20210726
                    dispoint1=0.    !20210726
                    nintf=backf(igapbf)%relat(kpoin)%nintf
                    do iintf=1,nintf
                        !jtotv=trans(itotv)%listf(iintf)  !20231026
                        jtotv=backf(igapbf)%relat(kpoin)%listf(iintf)    !20231026
                        dispoint1=dispoint1+result_zero(jtotv)*backf(igapbf)%relat(kpoin)%rintf(iintf)  !20231026
                    end do



                    qi=abs((backf(igapbf)%dism(kpoin,trstep)-dispoint1)/dispoint1)  !20220101
                    qerr=qerr+(backf(igapbf)%dism(kpoin,trstep)-dispoint1)**2 !20220101
                    qabs=qabs+dispoint1**2
                    write(7,20)kpoin,jdimn,qi,backf(igapbf)%dism(kpoin,trstep),dispoint1 !20220101

                end do
                tQerr=tQerr+Qerr;tQabs=tQabs+Qabs
                write(7,*)'step residual displacements (sum of squre),      Qerr=',qerr
                write(7,*)'step          displacements (sum of squre),      Qabs=',qabs
                write(7,*)'step  relative  displacementerror,root of (qerr/Qabs)=',sqrt(qerr/qabs)
            end do


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


        end do     !! loop for istep
        !if(kstab==0.and.(nforce/=0.or.ngaps/=0))call write_force_interface
        if(cwater/=0.and.delgroup>0)deallocate(coef_water)
    end do !!iincs
    write(7,*)'total residual displacements (sum of squre),      tQerr=',tqerr
    write(7,*)'total          displacements (sum of squre),      tQabs=',tqabs
    write(7,*)'total  relative  displacementerror,root of (tqerr/tQabs)=',sqrt(tqerr/tqabs)

20  format(2I10, 3e15.5)
    !if(winit==-1) call out_next_write
    if(winit==-1*iblks) call out_next_write !20231215YULI

    !if(ngaps/=0)deallocate(tofor0)  !!ctt2005

    END SUBROUTINE back_analysis  !20150925

    SUBROUTINE back_d_analysis  !20210820

    logical logx
    character(80)text
    integer(ink) itotv,ielem,irst,trstep0,ipoin,idofn,ij,idofix,ldofix,idelgroup,i0,ipairs,k
    real   (irk) xtime,time_begin,detal,ttime0,coef,xij,qi,qerr,qabs
    real   (irk) dispoint0,dispoint1,tQerr,tQabs,dispoint1g,dispoint0g  !20210726
    real   (irk) dis2e,diser,kmodu

    real   (irk),allocatable::rvectorm(:),value(:)
    integer(ink) iintf,nintf,iieq,igapbf,mdism,jpoinx   !!int2000

    integer(ink) i,iincs_i,iblks_i,istep_i,inode,jnode,ivalue,ivalue_point,bblks  !20200819
    integer(ink),pointer::listf(:)  !20200819
    real   (irk),pointer::rintf(:)  !20200819


    integer(ink) igapb,npgblock,jpoin,igaps,ipair,idimn,itotvbt,jdimn, &  !!ctt2005
        jtotv,kpoin,lpoin,jtotvbt,npairs,cwater,jgaps,jpair,kdimn,jpoin0,itotv0   !!ctt2005
    real   (irk),allocatable::rot(:,:),tofor0(:),uireact(:,:)  !!ctt2005
    real   (irk),allocatable::unitl(:),unitg(:),cmatrixl(:,:) !!ctt2005

    if (meshc==1.or.rmesh/=0)rewind(mainunit)
    if(Bparameter/=0.and.iblks==1)rewind(mainunit)  !20190810

    !read(observ_unit,*)text


    read(mainunit,*)text
    read(mainunit,*)nincs

    print *,' in back_d_analysis,trstep=',trstep
    tQerr=0.;tQabs=0.

    if(ngaps/=0.or.nrcsteel/=0)allocate(tofor0(ntotv)) !!ctt2005

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


        ttime0=ttime
        trstep0=trstep
        !if(nbackf/=0) & !20210726,20230430
        !read(mainunit,*)text

        do istep=inc_step,nstep,inc_step

            if(iblks>=stab_matde)call stab_initialize

            write(chkunit,*)'Increment step=',istep
            if(nbackf/=0)trstep=trstep0+istep  !20220101
            if(outintr>0.and.iblks>=outintr)trstep=trstep0+istep  !20200226
            xtime=ditime*istep
            ttime=ttime0+ditime*istep !! only for output

            result_zero=0.0  !201605
            do igapb=1,ngapb
                if(gapb(igapb)%nrdof/=0)gapb(igapb)%rdisp_zero=0.
            end do


            do igaps=1,ngaps   !20190810
                npairs=gaps(igaps)%npairs
                do ipairs=1,npairs
                    if(gaps(igaps)%pair_process(ipairs)==0)cycle  !20200331
                    gaps(igaps)%dxyz0(:,ipairs)=0.
                    gaps(igaps)%dxyz(:,ipairs)=0.
                    gaps(igaps)%ctforce0(:,ipairs)=0.
                    gaps(igaps)%ctforce(:,ipairs)=0.
                end do
            end do   !20190810
            call gpvar_initial  !201605

            if(nbackf/=0)then  !20210804
                do idofn=1,nbackf  !20210804
                    mdism=backf(idofn)%mdism
                    backf(idofn)%ic=1
                    do i0=1,mdism
                        if(abs(backf(idofn)%dism(i0,trstep)-999.)<.01) & !20220101
                            backf(idofn)%ic(i0)=0
                    end do
                enddo !20210804
            endif  !20210804


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
                if(type_load=='DISCONTROL')preact0=prescrib(1)%rdofix
                if((ngaps/=0.or.nrcsteel/=0).and.mdiv==1)tofor0=tofor !!ctt2005
                if((ngaps/=0.or.nrcsteel/=0).and.mdiv/=1)tofor0=toform !!ctt2005
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


1                   if(type_load=='LOAD2'.or.(kstat==2.and.iiter.le.2).or.(type_load/='LOAD2'.and.kstat/=2.and.iiter==1).or.  &
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

                        call eload_initialize

                        call residu_f

                        if(rmesh>0.and.nelem1>0)call residu_f1
                        if(rmesh>1.and.nelem2>0)call residu_f2
                        call eload_field
                        if(ground_inf/=0)call semi_inf_load
                        call force_internal
                    endif !for iiter==1 and istep==inc_step .and.idiv==1  temperature


                    if(ngaps/=0.and.iblks>=iblks_bt.and.iiter==1.and.mdiv==1)call ctfor_to_tofor(tofor0,tofor)  !!ctt2005
                    if(ngaps/=0.and.iblks>=iblks_bt.and.iiter==1.and.mdiv/=1)call ctfor_to_tofor(tofor0,toform)  !!ctt2005
                    if(nrcsteel/=0.and.iiter==1.and.mdiv==1)call csfor_to_tofor(tofor0,tofor)  !!20210328
                    if(nrcsteel/=0.and.iiter==1.and.mdiv/=1)call csfor_to_tofor(tofor0,toform)  !!20210328

                    if (mdiv/=1) then
                        if(idiv==1.and.iiter==1.and.allocated(torel))toform=toform+torel  !!20210328
                    else
                        if(iiter==1.and.allocated(torel))tofor=tofor+torel !!20210328
                    endif


                    if(neuman==1.and.(istep/=1.or.iiter/=1).and.kresl/=0)  goto 2  !ctt2005 , change position!
                    if(type_nl==8.and.(iiter>1.or.(kstat==2.and.iiter>2))) goto 2  !MNR
                    if ((type_solver=='PROFILE'.or.type_solver=='PARDISO').and.kresl/=0)then
                        operation='FACTORIZE'
                        call solve
                    end if
2                   continue

                    logx=ngaps/=0.and.(iiter==1.and.istep==inc_step.and.iincs==(lincs+1)).and.iblks==iblks_bt
                    if (logx)then !ctt2005
                        if (restart_ctt==0)then !restart_ctt
                            kdimn=ndimn
                            if(block_stab==1)kdimn=3*(ndimn-1)
                            allocate(rot(kdimn,kdimn))
                            rot=0.
                            do igapbf=1,nbackf
                                igapb=backf(igapbf)%groupb
                                npgblock=gapb(igapb)%npgblock
                                if(backf(igapbf)%mdism>npgblock*kdimn) allocate(uireact(backf(igapbf)%mdism,npgblock*kdimn))
                                gapb(igapb)%cmatrix=0.
                                do ipoin=1,npgblock
                                    igaps=gapb(igapb)%nodegblock_igaps(ipoin)
                                    ipair=gapb(igapb)%nodegblock_ipairs(ipoin)
                                    ij=gapb(igapb)%nodegblock_onetwo(ipoin)
                                    jpoinx=gaps(igaps)%pairnode(ij,ipair)
                                    coef=1.
                                    !if(ij==2)coef=-1.  !20210820
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
                                        unitg=unitl  !20210820
                                        !unitg=transpose(rot).x.unitl
                                        rvector=0.

                                        call  unit_dis_force_trans(igapb,jpoinx,kdimn,unitg,rvector)

                                        operation='SOLVE'
                                        call solve
                                        !write(7,*)'itotvbt=',itotvbt
                                        !write(7,*)'result=',result
                                        if(backf(igapbf)%mdism==npgblock*kdimn)then
                                            do kpoin=1,backf(igapbf)%mdism
                                                !jdimn=backf(igapbf)%listdim(kpoin)
                                                dispoint1=0.    !20230523
                                                nintf=backf(igapbf)%relat(kpoin)%nintf
                                                do iintf=1,nintf
                                                    jtotv=backf(igapbf)%relat(kpoin)%listf(iintf)
                                                    dispoint1=dispoint1+(result_zero(jtotv)+result(jtotv))*backf(igapbf)%relat(kpoin)%rintf(iintf)
                                                end do

                                                jtotvbt=kpoin
                                                !if(jpoin0==0)then
                                                gapb(igapb)%cmatrix(jtotvbt,itotvbt)=dispoint1
                                                !else
                                                !gapb(igapb)%cmatrix(jtotvbt,itotvbt)=dispoint1-dispoint0
                                                !endif
                                            end do  !kpoin
                                        else if(backf(igapbf)%mdism>npgblock*kdimn)then
                                            do kpoin=1,backf(igapbf)%mdism
                                                !jdimn=backf(igapbf)%listdim(kpoin)
                                                dispoint1=0.    !20230523
                                                nintf=backf(igapbf)%relat(kpoin)%nintf
                                                do iintf=1,nintf
                                                    jtotv=backf(igapbf)%relat(kpoin)%listf(iintf)
                                                    dispoint1=dispoint1+(result_zero(jtotv)+result(jtotv))*backf(igapbf)%relat(kpoin)%rintf(iintf)
                                                    !write(7,*)'kpoin=',kpoin,'jtotv=',jtotv,'result=',result(jtotv),'rintf=',backf(igapbf)%relat(kpoin)%rintf(iintf)
                                                end do

                                                jtotvbt=kpoin
                                                !if(jpoin0==0)then
                                                uireact(jtotvbt,itotvbt)=dispoint1  !20210726

                                                !write(7,*)'jtotvbt,itotvbt=',jtotvbt,itotvbt,'uireact=',uireact(jtotvbt,itotvbt)
                                                !   else
                                                !uireact(jtotvbt,itotvbt)=dispoint1-dispoint0    !20210726
                                                !endif
                                            end do  !kpoin

                                        endif
                                    end do  !idimn
                                    deallocate(unitl,unitg)
                                end do  !ipoin

                                if(backf(igapbf)%mdism>npgblock*kdimn)then
                                    gapb(igapb)%uireact=uireact
                                    gapb(igapb)%cmatrix(1:npgblock*kdimn,1:npgblock*kdimn)=transpose(uireact).x.gapb(igapb)%uireact
                                    deallocate(uireact)
                                endif


                            end do  !igapb
                            !call forAdirect_back_analysis !fzx !形成A矩阵


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
                            call forAdirect_back_analysis !fzx !形成A矩阵
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
                    rvector=0.0
                    !write(7,*)'iiter=',iiter,'itotv,ieq,tofor,stfor,rvector='
                    if (type_solver/='JPCG') then
                        do itotv=1,ntotv
                            if (totveq(itotv)/=0)then
                                if (mdiv/=1)then
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                        toform(itotv)-stfor(itotv)

                                else
                                    rvector(totveq(itotv))=rvector(totveq(itotv))+ &
                                        tofor(itotv)-stfor(itotv)
                                    !if(abs(rvector(totveq(itotv)))>1.e-3) &
                                    !                   write(7,*)itotv,totveq(itotv),tofor(itotv),stfor(itotv), rvector(totveq(itotv))
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


                    if (type_nl==8)then
                        if(kstat/=2)call bfgsr(iiter)
                        if(kstat==2)call bfgsr(iiter-1)
                    else
                        operation='SOLVE'
                        call solve
                    endif
                    result_zero_e=result_zero_e+result   !20210820


                    if(ngaps/=0.and.iblks>=abs(iblks_bt))call solve_back_d_analysis  !!ctt2005



                    call TIME(char_time)
                    print *, 'time: ', char_time
                    write(chkunit,*)'time: ', char_time


                    call varupdate
                    call eload_initialize
                    call residu_f
                    if(rmesh>0.and.nelem1>0)call residu_f1
                    if(rmesh>1.and.nelem2>0)call residu_f2
                    if(type_load/='LOAD2'.or.(type_load=='LOAD2'.and.idiv==2))then   !806
                        call eload_field
                        if(ground_inf/=0)call semi_inf_load

                        call reaction_prescribed


                        call conver_load
                        if(nchek==0) call conver_nodal_value

                        if(nchek==0)exit !tcl
                    endif !ep2010
10                  continue
                    print *,'miter=',miter,'iiter=',iiter
                end do   !! loop for iiter

                if(type_load/='LOAD2') &
                    call gpvarupdate

                if(rmesh>0.and.nelem1>0)call gpvarupdate1
                if(rmesh>1.and.nelem2>0)call gpvarupdate2
            end do    !! for idiv
            if(type_load=='LOAD2') &
                call gpvarupdate
            if(modf_dis_blocks(iblks)==1)call construction_dis_modify


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

            do igapbf=1,nbackf
                write(7,*)'nodal number,direction,ratio,foundation_dis,dam_dis'
                dis2e=0.;diser=0.
                qi=0.
                do kpoin=1,backf(igapbf)%mdism
                    if(backf(igapbf)%ic(kpoin)==0)cycle  !20210726
                    jdimn=backf(igapbf)%listdim(kpoin)
                    dispoint1=0.;dispoint1g=0.   !20230523  20231026
                    nintf=backf(igapbf)%relat(kpoin)%nintf
                    do iintf=1,nintf
                        jtotv=backf(igapbf)%relat(kpoin)%listf(iintf)
                        dispoint1=dispoint1+result_zero_e(jtotv)*backf(igapbf)%relat(kpoin)%rintf(iintf)
                        dispoint1g=dispoint1g+(result_zero(jtotv)-result_zero_e(jtotv))*backf(igapbf)%relat(kpoin)%rintf(iintf)
                    end do

                    dispoint1g=backf(igapbf)%dism(kpoin,trstep)-dispoint1  !20231026

                    qi=dispoint1g/dispoint1
                    !dis2e=dis2e+dispoint1**2
                    !diser=diser+(backf(igapbf)%dism(kpoin,trstep)-dispoint1g)*dispoint1

                    write(7,20)kpoin,jdimn,qi,dispoint1g,dispoint1
                end do
            end do

            !    kmodu=diser/dis2e
            !write(7,*)'kmodu=',kmodu

            do igapbf=1,nbackf
                write(7,*)'relative displacement error/nodal number,direction,relative error,u_observation,u_computation'
                qi=0.; Qerr=0.;Qabs=0.
                do kpoin=1,backf(igapbf)%mdism
                    if(backf(igapbf)%ic(kpoin)==0)cycle  !20210726
                    jdimn=backf(igapbf)%listdim(kpoin)
                    dispoint1=0.    !20230523
                    nintf=backf(igapbf)%relat(kpoin)%nintf
                    do iintf=1,nintf
                        jtotv=backf(igapbf)%relat(kpoin)%listf(iintf)
                        dispoint1=dispoint1+(result_zero(jtotv)+result(jtotv))*backf(igapbf)%relat(kpoin)%rintf(iintf)
                    end do


                    qi=abs((backf(igapbf)%dism(kpoin,trstep)-dispoint1)/dispoint1) !20220101
                    qerr=qerr+(backf(igapbf)%dism(kpoin,trstep)-dispoint1)**2  !20220101
                    qabs=qabs+dispoint1**2
                    write(7,20)kpoin,jdimn,qi,backf(igapbf)%dism(kpoin,trstep),dispoint1 !20220101

                end do
                tQerr=tQerr+Qerr;tQabs=tQabs+Qabs
                write(7,*)'step residual displacements (sum of squre),      Qerr=',qerr
                write(7,*)'step          displacements (sum of squre),      Qabs=',qabs
                write(7,*)'step  relative  displacementerror,root of (qerr/Qabs)=',sqrt(qerr/qabs)
            end do


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


        end do     !! loop for istep
        !if(kstab==0.and.(nforce/=0.or.ngaps/=0))call write_force_interface
        if(cwater/=0.and.delgroup>0)deallocate(coef_water)
    end do !!iincs
    write(7,*)'total residual displacements (sum of squre),      tQerr=',tqerr
    write(7,*)'total          displacements (sum of squre),      tQabs=',tqabs
    write(7,*)'total  relative  displacementerror,root of (tqerr/tQabs)=',sqrt(tqerr/tqabs)

20  format(2I10, 3e15.5)
    !if(winit==-1) call out_next_write
    if(winit==-1*iblks) call out_next_write !20231215YULI
    !if(ngaps/=0)deallocate(tofor0)  !!ctt2005

    END SUBROUTINE back_d_analysis  !20210820
