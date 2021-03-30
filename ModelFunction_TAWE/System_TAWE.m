function [TF,Vel,Cor,Jac,TransDis,Quat,MM,GFI,GFF,GFU,JacU,JacCons,CorCons,GFCons,CenMat,CorMatLeft,CorMatRight]=System_TAWE(t,q,p,u,s)
    
    [TF,Vel,Cor,Jac,TransDis,Quat]=numKinematics_TAWE(t,q,p,u,s);
    [MM,GFI,CenMat,CorMatLeft,CorMatRight]=Inertial_TAWE(t,q,p,u,s,TF,Vel,Cor,Jac);
    [GFF]=Force_TAWE(t,q,p,u,s,TF,TransDis,Vel,Jac,Quat);
    [GFU,JacU]=Input_TAWE(t,q,p,u,s,TF,TransDis,Vel,Jac,Quat);
    [JacCons,CorCons,GFCons]=Constraint_TAWE(t,q,p,u,s,TransDis,Vel,Cor,Jac,Quat);

	function [frameTF,frameVel,frameCorAcc,frameJacobian,frameTransDis,frameRotQuat]=numKinematics_TAWE(t,q,p,u,s)
	    LinkTF=linkTF(t,q,p,u,s);
	    LinkVel=linkVel(t,q,p,u,s);
	    LinkCorAcc=linkCorAcc(t,q,p,u,s);
	    LinkJacobian=linkJacobian(t,q,p,u,s);
	    LinkRotQuat=linkRotQuat(t,q,p,u,s);
	    FrameNum=26;
	    DOF=numel(q)/2;
	    
	    LinkTF=reshape((LinkTF),[4,4,FrameNum]);
	    LinkJacobian=reshape((LinkJacobian),[6,DOF,FrameNum]);
	    
	    frameTF=zeros(4,4,FrameNum);
	    frameTransDis=zeros(3,FrameNum);
	    frameVel=zeros(6,FrameNum);
	    frameCorAcc=zeros(6,FrameNum);
	    frameJacobian=zeros(6,DOF,FrameNum);
	    frameRotQuat=zeros(4,FrameNum);
	    
	    frameTF(:,:,1)=LinkTF(:,:,1);
	    frameTransDis(:,1)=frameTF(1:3,4,1);
	    frameRotQuat(:,1)=LinkRotQuat(:,1);
	    frameVel(:,1)=LinkVel(:,1);
	    frameCorAcc(:,1)=LinkCorAcc(:,1);
	    frameJacobian(:,:,1)=LinkJacobian(:,:,1);
	    frameCheckList=zeros(1,FrameNum);
	    frameCheckList(1)=1;
	    for frame=1:FrameNum
	        framePath=0;
	%SWITCHCASE_
			switch frame
			    case 1
			        framePath=[1];
			    case 2
			        framePath=[1  2];
			    case 3
			        framePath=[1  2  3];
			    case 4
			        framePath=[1  2  3  4];
			    case 5
			        framePath=[1  2  3  4  5];
			    case 6
			        framePath=[1  2  3  4  5  6];
			    case 7
			        framePath=[1  2  3  4  5  6  7];
			    case 8
			        framePath=[1  2  3  4  5  6  8];
			    case 9
			        framePath=[1  2  3  4  5  6  9];
			    case 10
			        framePath=[1   2   3   4   5   6  10];
			    case 11
			        framePath=[1   2   3   4   5   6  11];
			    case 12
			        framePath=[1   2  12];
			    case 13
			        framePath=[1   2  12  13];
			    case 14
			        framePath=[1   2  12  13  14];
			    case 15
			        framePath=[1   2  12  13  14  15];
			    case 16
			        framePath=[1   2  12  13  14  15  16];
			    case 17
			        framePath=[1   2  12  13  14  15  16  17];
			    case 18
			        framePath=[1   2  12  13  14  15  16  17  18];
			    case 19
			        framePath=[1   2  12  13  14  15  16  17  18  19];
			    case 20
			        framePath=[1   2  12  13  14  15  16  17  18  19  20];
			    case 21
			        framePath=[1   2  12  13  14  21];
			    case 22
			        framePath=[1   2  12  13  14  15  22];
			    case 23
			        framePath=[1   2  12  13  14  15  16  23];
			    case 24
			        framePath=[1   2  12  13  14  15  16  17  24];
			    case 25
			        framePath=[1   2  12  13  14  15  16  17  18  25];
			    case 26
			        framePath=[1   2  12  13  14  15  16  17  18  19  26];
			end
	        [frameCheckList,LinkTF,frameTF,LinkVel,frameVel,LinkCorAcc,frameCorAcc,LinkJacobian,frameJacobian,frameTransDis,LinkRotQuat,frameRotQuat]=frameKineRecur(framePath,frameCheckList,LinkTF,frameTF,LinkVel,frameVel,LinkCorAcc,frameCorAcc,LinkJacobian,frameJacobian,frameTransDis,LinkRotQuat,frameRotQuat);
	    end
	    function output_=skew3(w_)
	        output_=[0 -w_(3) w_(2) ; w_(3) 0 -w_(1) ; -w_(2) w_(1) 0 ];
	    end
	    function p_ = quatMultiply(q_, r_)
	        p_ = [q_(1).*r_(1) - q_(2).*r_(2) - q_(3).*r_(3) - q_(4).*r_(4);
	             q_(1).*r_(2) + r_(1).*q_(2) + q_(3).*r_(4) - q_(4).*r_(3);
	             q_(1).*r_(3) + r_(1).*q_(3) + q_(4).*r_(2) - q_(2).*r_(4);
	             q_(1).*r_(4) + r_(1).*q_(4) + q_(2).*r_(3) - q_(3).*r_(2)];
	    end
	    function [frameCheckList_,TF_,frameTF_,Vel_,frameVel_,CorAcc_,frameCorAcc_,Jacobian_,frameJacobian_,frameTransDis_,RotQuat_,frameRotQuat_]=frameKineRecur(framePath_,frameCheckList_,TF_,frameTF_,Vel_,frameVel_,CorAcc_,frameCorAcc_,Jacobian_,frameJacobian_,frameTransDis_,RotQuat_,frameRotQuat_)
	        curNum_=framePath_(end);
	        if(~frameCheckList_(curNum_))
	            preNum_=framePath_(end-1);
	            if(frameCheckList_(preNum_))
	                frameTF_(:,:,curNum_)=frameTF_(:,:,preNum_)*TF_(:,:,curNum_);
	                frameTransDis_(:,curNum_)=frameTF_(1:3,4,curNum_);
	                frameRotQuat_(:,curNum_)=quatMultiply(frameRotQuat_(:,preNum_),RotQuat_(:,curNum_));
	                frameRotQuat_(:,curNum_)=frameRotQuat_(:,curNum_)./norm(frameRotQuat_(:,curNum_));
	                curTransVel=Vel_(1:3,curNum_);
	                curAngVel=Vel_(4:6,curNum_);
	                curCorTransAcc=CorAcc_(1:3,curNum_);
	                curCorAngAcc=CorAcc_(4:6,curNum_);
	                curTransJacobian=Jacobian_(1:3,:,curNum_);
	                curAngJacobian=Jacobian_(4:6,:,curNum_);
	                preTransVel=frameVel_(1:3,preNum_);
	                preAngVel=frameVel_(4:6,preNum_);
	                preCorTransAcc=frameCorAcc_(1:3,preNum_);
	                preCorAngAcc=frameCorAcc_(4:6,preNum_);
	                preTransJacobian=frameJacobian_(1:3,:,preNum_);
	                preAngJacobian=frameJacobian_(4:6,:,preNum_);
	                
	                preRotMat=frameTF_(1:3,1:3,preNum_);
	                preRotCurTransDis=preRotMat*TF_(1:3,4,curNum_);
	                
	                frameVel_(1:3,curNum_)=preTransVel+preRotMat*curTransVel+cross(preAngVel,preRotCurTransDis);
	                frameVel_(4:6,curNum_)=preAngVel+preRotMat*curAngVel;
	                frameCorAcc_(1:3,curNum_)=preCorTransAcc+cross(preAngVel,cross(preAngVel,preRotCurTransDis))...
	                                        +preRotMat*curCorTransAcc+cross(preCorAngAcc,preRotCurTransDis)...
	                                        +2*cross(preAngVel,preRotMat*curTransVel);
	                frameCorAcc_(4:6,curNum_)=preCorAngAcc+preRotMat*curCorAngAcc+cross(preAngVel,preRotMat*curAngVel);
	                %skewTransDis=[0 -preRotCurTransDis(3) preRotCurTransDis(2) ; preRotCurTransDis(3) 0 -preRotCurTransDis(1) ; -preRotCurTransDis(2) preRotCurTransDis(1) 0 ];
	                skewTransDis=skew3(preRotCurTransDis);
	                frameJacobian_(1:3,:,curNum_)=preTransJacobian+preRotMat*curTransJacobian-skewTransDis*preAngJacobian;
	                frameJacobian_(4:6,:,curNum_)=preAngJacobian+preRotMat*curAngJacobian;
	                frameCheckList_(curNum_)=1;
	            else
	                [frameCheckList_,TF_,frameTF_,Vel_,frameVel_,CorAcc_,frameCorAcc_,Jacobian_,frameJacobian_,frameTransDis_,RotQuat_,frameRotQuat_]=frameKineRecur(framePath_(1:end-1),frameCheckList_,TF_,frameTF_,Vel_,frameVel_,CorAcc_,frameCorAcc_,Jacobian_,frameJacobian_,frameTransDis_,RotQuat_,frameRotQuat_);
	            end
	        end
	    end
		function out1 = linkTF(t,in2,in3,in4,in5)
		%LINKTF
		%    OUT1 = LINKTF(T,IN2,IN3,IN4,IN5)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:15:17
		l21 = in3(19,:);
		l22 = in3(20,:);
		l23 = in3(21,:);
		l24 = in3(22,:);
		l25 = in3(23,:);
		l11x = in3(7,:);
		l11y = in3(8,:);
		l12x = in3(10,:);
		l11z = in3(9,:);
		l12y = in3(11,:);
		l13x = in3(13,:);
		l12z = in3(12,:);
		l13y = in3(14,:);
		l13z = in3(15,:);
		l21x = in3(16,:);
		l21y = in3(17,:);
		l21z = in3(18,:);
		lm13x = in3(42,:);
		lm13y = in3(43,:);
		lm13z = in3(44,:);
		lm22x = in3(45,:);
		lm22y = in3(46,:);
		lm23x = in3(48,:);
		lm22z = in3(47,:);
		lm23y = in3(49,:);
		lm24x = in3(51,:);
		lm23z = in3(50,:);
		lm24y = in3(52,:);
		lm25x = in3(54,:);
		lm24z = in3(53,:);
		lm25y = in3(55,:);
		lm26x = in3(57,:);
		lm25z = in3(56,:);
		lm26y = in3(58,:);
		lm27x = in3(60,:);
		lm26z = in3(59,:);
		lm27y = in3(61,:);
		lm27z = in3(62,:);
		q1dev__dt_0_ = in2(1,:);
		q1flex__dt_0_ = in2(2,:);
		q1sup__dt_0_ = in2(3,:);
		q2w1__dt_0_ = in2(4,:);
		q2w2__dt_0_ = in2(5,:);
		q2w3__dt_0_ = in2(6,:);
		q2w4__dt_0_ = in2(7,:);
		q2w5__dt_0_ = in2(8,:);
		q2w6__dt_0_ = in2(9,:);
		rTetra = in3(6,:);
		rf1IMU1x = in3(24,:);
		rf1IMU1y = in3(25,:);
		rf1IMU2x = in3(30,:);
		rf1IMU1z = in3(26,:);
		rf1IMU2y = in3(31,:);
		rf1IMU2z = in3(32,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf2Medz = in3(35,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		rf1Devz = in3(29,:);
		riw1 = in3(36,:);
		riw2 = in3(37,:);
		riw3 = in3(38,:);
		riw4 = in3(39,:);
		riw5 = in3(40,:);
		riw6 = in3(41,:);
		t2 = (cos(rf1IMU1z));
		t3 = (sin(rf1IMU1x));
		t4 = (sin(rf1IMU1z));
		t5 = (cos(rf1IMU1x));
		t6 = (sin(rf1IMU1y));
		t7 = (cos(rf1Devz));
		t8 = (sin(rf1Devx));
		t9 = (sin(rf1Devz));
		t10 = (cos(rf1Devx));
		t11 = (sin(rf1Devy));
		t12 = (cos(q1dev__dt_0_));
		t13 = (sin(q1dev__dt_0_));
		t14 = (sin(q1flex__dt_0_));
		t15 = (cos(q1flex__dt_0_));
		t16 = (sin(q1sup__dt_0_));
		t17 = (cos(rf1IMU2z));
		t18 = (sin(rf1IMU2x));
		t19 = (sin(rf1IMU2z));
		t20 = (cos(rf1IMU2x));
		t21 = (sin(rf1IMU2y));
		t22 = (cos(rf2Medz));
		t23 = (sin(rf2Medx));
		t24 = (sin(rf2Medz));
		t25 = (cos(rf2Medx));
		t26 = (sin(rf2Medy));
		t27 = (q2w1__dt_0_+riw1);
		t28 = (q2w5__dt_0_+riw5);
		t29 = (q2w6__dt_0_+riw6);
		t30 = (cos(rf1IMU2y));
		t31 = (t17.*t30);
		t32 = (cos(rf1IMU1y));
		t35 = (rf1IMU1x./2.0);
		t33 = (cos(t35));
		t37 = (rf1IMU1z./2.0);
		t34 = (cos(t37));
		t36 = (t33.^2);
		t38 = (t34.^2);
		t39 = (rf1IMU1y./2.0);
		t40 = (cos(rf1Devy));
		t43 = (rf1Devx./2.0);
		t41 = (cos(t43));
		t45 = (rf1Devz./2.0);
		t42 = (cos(t45));
		t44 = (t41.^2);
		t46 = (t42.^2);
		t47 = (rf1Devy./2.0);
		t48 = (cos(q1sup__dt_0_));
		t51 = (q1dev__dt_0_./2.0);
		t49 = (cos(t51));
		t53 = (q1flex__dt_0_./2.0);
		t50 = (cos(t53));
		t52 = (t49.^2);
		t54 = (t50.^2);
		t55 = (q1sup__dt_0_./2.0);
		t56 = (t19.*t30);
		t59 = (rf1IMU2x./2.0);
		t57 = (cos(t59));
		t61 = (rf1IMU2z./2.0);
		t58 = (cos(t61));
		t60 = (t57.^2);
		t62 = (t58.^2);
		t63 = (rf1IMU2y./2.0);
		t64 = (sqrt(6.0));
		t65 = (cos(rf2Medy));
		t68 = (rf2Medx./2.0);
		t66 = (cos(t68));
		t70 = (rf2Medz./2.0);
		t67 = (cos(t70));
		t69 = (t66.^2);
		t71 = (t67.^2);
		t72 = (rf2Medy./2.0);
		t73 = (sin(t27));
		t74 = (cos(t27));
		t75 = (q2w2__dt_0_+riw2);
		t76 = (q2w3__dt_0_+riw3);
		t77 = (q2w4__dt_0_+riw4);
		t78 = (sin(t29));
		t79 = (cos(t29));
		t80 = (t17.*t18.*t21);
		t81 = (t80-t19.*t20);
		t82 = (t60.*t62.*4.0);
		t83 = (cos(t63));
		t84 = (sin(t59));
		t85 = (sin(t63));
		t86 = (sin(t61));
		t87 = (t57.*t58.*t83.*t84.*t85.*t86.*8.0);
		t88 = (t60.*-2.0-t62.*2.0+t82+t87+1.0);
		t89 = (t18.*t30);
		t90 = (sqrt(3.0));
		t91 = (sin(t75));
		t92 = (cos(t75));
		t93 = (sin(t76));
		t94 = (cos(t76));
		t95 = (sin(t77));
		t96 = (cos(t77));
		t97 = (sin(t28));
		t98 = (cos(t28));
		t99 = (t18.*t19);
		t100 = (t17.*t20.*t21);
		t101 = (t99+t100);
		t102 = (t19.*t20.*t21);
		t103 = (t102-t17.*t18);
		t104 = (t20.*t30);
		outSize=[4,104];
		elementRow=[1,2,3,4,1,2,3,1,2,3,1,2,3,4,1,2,3,1,2,3,1,2,3,1,2,3,4,1,2,3,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,2,3,4,1,2,3,2,4,1,2,3,4,1,2,3,1,2,3,1,2,3,1,2,3,4,1,2,1,2,3,4,1,2,3,2,3,4,1,2,3,2,3,2,4,1,2,3,2,3,2,4,1,3,2,1,3,2,4,1,2,1,2,3,1,3,4,1,2,3,1,2,3,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4];
		elementCol=[1,2,3,4,5,5,5,6,6,6,7,7,7,8,9,9,9,10,10,10,11,11,11,12,12,12,12,13,13,13,14,14,14,15,15,15,16,17,18,19,20,20,20,20,21,21,21,22,22,22,23,23,23,24,25,26,27,28,28,28,28,29,30,31,32,32,32,32,33,34,35,36,36,36,36,37,38,39,40,40,40,41,42,43,44,44,45,46,47,48,49,49,49,50,50,50,51,51,51,52,52,52,52,53,53,54,54,55,56,57,58,58,59,59,60,61,62,62,63,63,64,64,65,66,66,67,67,68,68,69,69,70,71,71,72,72,73,73,74,74,75,76,76,76,77,77,77,78,78,78,79,79,79,80,80,80,80,81,82,83,84,84,84,84,85,86,87,88,88,88,88,89,90,91,92,92,92,92,93,94,95,96,96,96,96,97,98,99,100,100,100,100,101,102,103,104,104,104,104];
		elementList=[1.0,1.0,1.0,1.0,t2.*t32,t4.*t32,-t6,-t4.*t5+t2.*t3.*t6,t36.*-2.0-t38.*2.0+t36.*t38.*4.0+t33.*t34.*cos(t39).*sin(t35).*sin(t37).*sin(t39).*8.0+1.0,t3.*t32,t3.*t4+t2.*t5.*t6,-t2.*t3+t4.*t5.*t6,t5.*t32,1.0,t7.*t40,t9.*t40,-t11,-t9.*t10+t7.*t8.*t11,t44.*-2.0-t46.*2.0+t44.*t46.*4.0+t41.*t42.*cos(t47).*sin(t43).*sin(t45).*sin(t47).*8.0+1.0,t8.*t40,t8.*t9+t7.*t10.*t11,-t7.*t8+t9.*t10.*t11,t10.*t40,l11x,l11y,l11z,1.0,t12.*t48,t13.*t48,-t16,-t13.*t15+t12.*t14.*t16,t52.*-2.0-t54.*2.0+t52.*t54.*4.0+t49.*t50.*cos(t55).*sin(t51).*sin(t53).*sin(t55).*8.0+1.0,t14.*t48,t13.*t14+t12.*t15.*t16,-t12.*t14+t13.*t15.*t16,t15.*t48,1.0,1.0,1.0,1.0,l12x,l12y,l12z,1.0,t31,t56,-t21,t81,t88,t89,t101,t103,t104,1.0,1.0,1.0,1.0,lm13x,lm13y,lm13z,1.0,1.0,1.0,1.0,-rTetra,rTetra.*t64.*(-1.0./6.0),rTetra.*t90.*(-1.0./3.0),1.0,1.0,1.0,1.0,rTetra,rTetra.*t64.*(-1.0./6.0),rTetra.*t90.*(-1.0./3.0),1.0,1.0,1.0,1.0,rTetra.*t64.*(-1.0./6.0),rTetra.*t90.*(2.0./3.0),1.0,1.0,1.0,1.0,(rTetra.*t64)./2.0,1.0,1.0,1.0,1.0,1.0,t22.*t65,t24.*t65,-t26,-t24.*t25+t22.*t23.*t26,t69.*-2.0-t71.*2.0+t69.*t71.*4.0+t66.*t67.*cos(t72).*sin(t68).*sin(t70).*sin(t72).*8.0+1.0,t23.*t65,t23.*t24+t22.*t25.*t26,-t22.*t23+t24.*t25.*t26,t25.*t65,l21x,l21y,l21z,1.0,t74,t73,-t73,t74,1.0,1.0,1.0,t92,t91,-t91,t92,1.0,1.0,t94,t93,-t93,t94,l21,1.0,1.0,t96,t95,-t95,t96,l22,1.0,t98,-t97,1.0,t97,t98,l23,1.0,t79,t78,-t78,t79,1.0,l24,l25,1.0,t31,t81,t101,t56,t88,t103,-t21,t89,t104,l13x,l13y,l13z,1.0,1.0,1.0,1.0,lm22x,lm22y,lm22z,1.0,1.0,1.0,1.0,lm23x,lm23y,lm23z,1.0,1.0,1.0,1.0,lm24x,lm24y,lm24z,1.0,1.0,1.0,1.0,lm25x,lm25y,lm25z,1.0,1.0,1.0,1.0,lm26x,lm26y,lm26z,1.0,1.0,1.0,1.0,lm27x,lm27y,lm27z,1.0];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = linkVel(t,in2,in3,in4,in5)
		%LINKVEL
		%    OUT1 = LINKVEL(T,IN2,IN3,IN4,IN5)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:15:14
		q1dev__dt_0_ = in2(1,:);
		q1dev__dt_1_ = in2(10,:);
		q1flex__dt_1_ = in2(11,:);
		q1sup__dt_0_ = in2(3,:);
		q1sup__dt_1_ = in2(12,:);
		q2w1__dt_1_ = in2(13,:);
		q2w2__dt_1_ = in2(14,:);
		q2w3__dt_1_ = in2(15,:);
		q2w4__dt_1_ = in2(16,:);
		q2w5__dt_1_ = in2(17,:);
		q2w6__dt_1_ = in2(18,:);
		t2 = cos(q1dev__dt_0_);
		t3 = cos(q1sup__dt_0_);
		t4 = sin(q1dev__dt_0_);
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-q1sup__dt_1_.*t4+q1flex__dt_1_.*t2.*t3,q1sup__dt_1_.*t2+q1flex__dt_1_.*t3.*t4,q1dev__dt_1_-q1flex__dt_1_.*sin(q1sup__dt_0_),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,q2w1__dt_1_,0.0,0.0,0.0,q2w2__dt_1_,0.0,0.0,0.0,0.0,0.0,q2w3__dt_1_,0.0,0.0,0.0,0.0,0.0,q2w4__dt_1_,0.0,0.0,0.0,0.0,0.0,0.0,q2w5__dt_1_,0.0,0.0,0.0,0.0,0.0,0.0,q2w6__dt_1_,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,26]);
		end
		function out1 = linkRotQuat(t,in2,in3,in4,in5)
		%LINKROTQUAT
		%    OUT1 = LINKROTQUAT(T,IN2,IN3,IN4,IN5)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:15:15
		q1dev__dt_0_ = in2(1,:);
		q1flex__dt_0_ = in2(2,:);
		q1sup__dt_0_ = in2(3,:);
		q2w1__dt_0_ = in2(4,:);
		q2w2__dt_0_ = in2(5,:);
		q2w3__dt_0_ = in2(6,:);
		q2w4__dt_0_ = in2(7,:);
		q2w5__dt_0_ = in2(8,:);
		q2w6__dt_0_ = in2(9,:);
		rf1IMU1x = in3(24,:);
		rf1IMU1y = in3(25,:);
		rf1IMU2x = in3(30,:);
		rf1IMU1z = in3(26,:);
		rf1IMU2y = in3(31,:);
		rf1IMU2z = in3(32,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf2Medz = in3(35,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		rf1Devz = in3(29,:);
		riw1 = in3(36,:);
		riw2 = in3(37,:);
		riw3 = in3(38,:);
		riw4 = in3(39,:);
		riw5 = in3(40,:);
		riw6 = in3(41,:);
		t2 = rf1IMU1x./2.0;
		t3 = rf1IMU1y./2.0;
		t4 = rf1IMU1z./2.0;
		t5 = rf1Devx./2.0;
		t6 = rf1Devy./2.0;
		t7 = rf1Devz./2.0;
		t8 = q1dev__dt_0_./2.0;
		t9 = q1flex__dt_0_./2.0;
		t10 = q1sup__dt_0_./2.0;
		t11 = rf1IMU2x./2.0;
		t12 = rf1IMU2y./2.0;
		t13 = rf1IMU2z./2.0;
		t14 = rf2Medx./2.0;
		t15 = rf2Medy./2.0;
		t16 = rf2Medz./2.0;
		t17 = sin(t11);
		t18 = sin(t12);
		t19 = sin(t13);
		t20 = t17.*t18.*t19;
		t21 = cos(t11);
		t22 = cos(t12);
		t23 = cos(t13);
		t24 = t21.*t22.*t23;
		t25 = t20+t24;
		t26 = cos(t3);
		t27 = cos(t4);
		t28 = sin(t2);
		t29 = cos(t2);
		t30 = sin(t3);
		t31 = sin(t4);
		t32 = cos(t6);
		t33 = cos(t7);
		t34 = sin(t5);
		t35 = cos(t5);
		t36 = sin(t6);
		t37 = sin(t7);
		t38 = cos(t8);
		t39 = cos(t10);
		t40 = sin(t9);
		t41 = cos(t9);
		t42 = sin(t8);
		t43 = sin(t10);
		t44 = cos(t15);
		t45 = cos(t16);
		t46 = sin(t14);
		t47 = cos(t14);
		t48 = sin(t15);
		t49 = sin(t16);
		t50 = q2w2__dt_0_./2.0;
		t51 = riw2./2.0;
		t52 = t50+t51;
		t53 = q2w3__dt_0_./2.0;
		t54 = riw3./2.0;
		t55 = t53+t54;
		t56 = q2w4__dt_0_./2.0;
		t57 = riw4./2.0;
		t58 = t56+t57;
		t59 = t17.*t22.*t23;
		t60 = q2w5__dt_0_./2.0;
		t61 = riw5./2.0;
		t62 = t60+t61;
		t63 = t18.*t21.*t23;
		t64 = t17.*t19.*t22;
		t65 = q2w1__dt_0_./2.0;
		t66 = riw1./2.0;
		t67 = t65+t66;
		t68 = q2w6__dt_0_./2.0;
		t69 = riw6./2.0;
		t70 = t68+t69;
		t71 = t19.*t21.*t22;
		out1 = reshape([1.0,0.0,0.0,0.0,t26.*t27.*t29+t28.*t30.*t31,t26.*t27.*t28-t29.*t30.*t31,t26.*t28.*t31+t27.*t29.*t30,-t27.*t28.*t30+t26.*t29.*t31,t32.*t33.*t35+t34.*t36.*t37,t32.*t33.*t34-t35.*t36.*t37,t32.*t34.*t37+t33.*t35.*t36,-t33.*t34.*t36+t32.*t35.*t37,t38.*t39.*t41+t40.*t42.*t43,t38.*t39.*t40-t41.*t42.*t43,t39.*t40.*t42+t38.*t41.*t43,-t38.*t40.*t43+t39.*t41.*t42,1.0,0.0,0.0,0.0,t25,t59-t18.*t19.*t21,t63+t64,t71-t17.*t18.*t23,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,t44.*t45.*t47+t46.*t48.*t49,t44.*t45.*t46-t47.*t48.*t49,t44.*t46.*t49+t45.*t47.*t48,-t45.*t46.*t48+t44.*t47.*t49,cos(t67),0.0,0.0,sin(t67),cos(t52),sin(t52),0.0,0.0,cos(t55),sin(t55),0.0,0.0,cos(t58),sin(t58),0.0,0.0,cos(t62),0.0,sin(t62),0.0,cos(t70),0.0,0.0,sin(t70),t25,-t59+t18.*t19.*t21,-t63-t64,-t71+t17.*t18.*t23,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0],[4,26]);
		end
		function out1 = linkCorAcc(t,in2,in3,in4,in5)
		%LINKCORACC
		%    OUT1 = LINKCORACC(T,IN2,IN3,IN4,IN5)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:15:15
		q1dev__dt_0_ = in2(1,:);
		q1dev__dt_1_ = in2(10,:);
		q1flex__dt_1_ = in2(11,:);
		q1sup__dt_0_ = in2(3,:);
		q1sup__dt_1_ = in2(12,:);
		t2 = cos(q1dev__dt_0_);
		t3 = sin(q1dev__dt_0_);
		t4 = sin(q1sup__dt_0_);
		t5 = cos(q1sup__dt_0_);
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-q1dev__dt_1_.*q1sup__dt_1_.*t2-q1dev__dt_1_.*q1flex__dt_1_.*t3.*t5-q1flex__dt_1_.*q1sup__dt_1_.*t2.*t4,-q1dev__dt_1_.*q1sup__dt_1_.*t3+q1dev__dt_1_.*q1flex__dt_1_.*t2.*t5-q1flex__dt_1_.*q1sup__dt_1_.*t3.*t4,-q1flex__dt_1_.*q1sup__dt_1_.*t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,26]);
		end
		function out1 = linkJacobian(t,in2,in3,in4,in5)
		%LINKJACOBIAN
		%    OUT1 = LINKJACOBIAN(T,IN2,IN3,IN4,IN5)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:15:18
		q1dev__dt_0_ = in2(1,:);
		q1sup__dt_0_ = in2(3,:);
		t2 = (cos(q1sup__dt_0_));
		t3 = (sin(q1dev__dt_0_));
		t4 = (cos(q1dev__dt_0_));
		outSize=[6,234];
		elementRow=[6,4,5,6,4,5,6,4,4,4,5,6];
		elementCol=[28,29,29,29,30,30,121,131,141,151,161,171];
		elementList=[1.0,t2.*t4,t2.*t3,-sin(q1sup__dt_0_),-t3,t4,1.0,1.0,1.0,1.0,1.0,1.0];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
	end

	function [inertM,inertGF,inertCenMat,inertCorMatLeft,inertCorMatRight]=Inertial_TAWE(t,q,p,u,s,TF_Global,Vel_Global,Cor_Global,Jac_Global)
	    BaseFrameList=[6   7  21  22  23  24  25  26];
	    MassList=Mass_TAWE(t,q,p,s,u);
	    MomentList=reshape(Moment_TAWE(t,q,p,s,u),[3 3 numel(MassList)]);
	    
	    inertCorMatLeft=zeros(numel(q)/2,6*numel(MassList));
	    inertCorMatRight=zeros(6*numel(MassList),numel(q)/2);
	    inertM=zeros(numel(q)/2,numel(q)/2,numel(MassList));
	    inertGF=zeros(numel(q)/2,numel(MassList));
	    inertCenMat=zeros(numel(q)/2,numel(q)/2,numel(MassList));
	    for bodyNum=1:numel(MassList)
	        rotor=TF_Global(1:3,1:3,BaseFrameList(bodyNum));
	        m=MassList(bodyNum);
	        I=rotor*MomentList(:,:,bodyNum)*rotor.';
	        w=Vel_Global(4:6,BaseFrameList(bodyNum));
	        a=Cor_Global(1:3,BaseFrameList(bodyNum));
	        alpha=Cor_Global(4:6,BaseFrameList(bodyNum));
	        vJacobian=Jac_Global(1:3,:,BaseFrameList(bodyNum));
	        wJacobian=Jac_Global(4:6,:,BaseFrameList(bodyNum));
	        
	        thisInertCorMatLeft=[m*vJacobian;I*wJacobian].';
	        thisInertCorMatRight=[vJacobian;wJacobian]; %Take the time derivative of this to calculate coriolis numerically
	        inertCorMatLeft(:,6*bodyNum+(-5:0))=thisInertCorMatLeft;
	        inertCorMatRight(6*bodyNum+(-5:0),:)=thisInertCorMatRight;
	        inertCenMat(:,:,bodyNum)=- wJacobian.'* skew3_InertDynamicsOnly_(I*w) * wJacobian; % On the same side of Mqddot
	        inertM(:,:,bodyNum)=thisInertCorMatLeft*thisInertCorMatRight;
	        inertGF(:,bodyNum)= - thisInertCorMatLeft*[a;alpha] - wJacobian.'* cross(w,I*w); % on the different side of Mqddot
	        % inertM(:,:,bodyNum)=vJacobian.'*m*vJacobian+wJacobian.'*I*wJacobian;
	        % inertGF(:,bodyNum)=-vJacobian.'*m*a-wJacobian.'*(I*alpha+cross(w,I*w)); % on the different side of Mqddot
	    end
	    function output_=skew3_InertDynamicsOnly_(w_)
	        output_=[0 -w_(3) w_(2) ; w_(3) 0 -w_(1) ; -w_(2) w_(1) 0 ];
	    end
		function out1 = Mass_TAWE(t,in2,in3,in4,in5)
		%MASS_TAWE
		%    OUT1 = MASS_TAWE(T,IN2,IN3,IN4,IN5)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:49
		m12 = in3(112,:);
		m13 = in3(63,:);
		m22 = in3(70,:);
		m23 = in3(77,:);
		m24 = in3(84,:);
		m25 = in3(91,:);
		m26 = in3(98,:);
		m27 = in3(105,:);
		out1 = [m12,m13,m22,m23,m24,m25,m26,m27];
		end
		function out1 = Moment_TAWE(t,in2,in3,in4,in5)
		%MOMENT_TAWE
		%    OUT1 = MOMENT_TAWE(T,IN2,IN3,IN4,IN5)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:50
		i12xx = in3(113,:);
		i12xy = in3(114,:);
		i13xx = in3(64,:);
		i12xz = in3(115,:);
		i12yy = in3(116,:);
		i13xy = in3(65,:);
		i12yz = in3(117,:);
		i13xz = in3(66,:);
		i13yy = in3(67,:);
		i12zz = in3(118,:);
		i13yz = in3(68,:);
		i13zz = in3(69,:);
		i22xx = in3(71,:);
		i22xy = in3(72,:);
		i23xx = in3(78,:);
		i22xz = in3(73,:);
		i22yy = in3(74,:);
		i23xy = in3(79,:);
		i24xx = in3(85,:);
		i22yz = in3(75,:);
		i23xz = in3(80,:);
		i23yy = in3(81,:);
		i24xy = in3(86,:);
		i25xx = in3(92,:);
		i22zz = in3(76,:);
		i23yz = in3(82,:);
		i24xz = in3(87,:);
		i24yy = in3(88,:);
		i25xy = in3(93,:);
		i26xx = in3(99,:);
		i23zz = in3(83,:);
		i24yz = in3(89,:);
		i25xz = in3(94,:);
		i25yy = in3(95,:);
		i26xy = in3(100,:);
		i27xx = in3(106,:);
		i24zz = in3(90,:);
		i25yz = in3(96,:);
		i26xz = in3(101,:);
		i26yy = in3(102,:);
		i27xy = in3(107,:);
		i25zz = in3(97,:);
		i26yz = in3(103,:);
		i27xz = in3(108,:);
		i27yy = in3(109,:);
		i26zz = in3(104,:);
		i27yz = in3(110,:);
		i27zz = in3(111,:);
		outSize=[3,24];
		elementRow=[1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3];
		elementCol=[1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,15,15,16,16,16,17,17,17,18,18,18,19,19,19,20,20,20,21,21,21,22,22,22,23,23,23,24,24,24];
		elementList=[i12xx,i12xy,i12xz,i12xy,i12yy,i12yz,i12xz,i12yz,i12zz,i13xx,i13xy,i13xz,i13xy,i13yy,i13yz,i13xz,i13yz,i13zz,i22xx,i22xy,i22xz,i22xy,i22yy,i22yz,i22xz,i22yz,i22zz,i23xx,i23xy,i23xz,i23xy,i23yy,i23yz,i23xz,i23yz,i23zz,i24xx,i24xy,i24xz,i24xy,i24yy,i24yz,i24xz,i24yz,i24zz,i25xx,i25xy,i25xz,i25xy,i25yy,i25yz,i25xz,i25yz,i25zz,i26xx,i26xy,i26xz,i26xy,i26yy,i26yz,i26xz,i26yz,i26zz,i27xx,i27xy,i27xz,i27xy,i27yy,i27yz,i27xz,i27yz,i27zz];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
	end

	function [CollectGF]=Force_TAWE(t,q,p,u,s,TF_Global,TransDis_Global,Vel_Global,Jac_Global,Quat_Global)
	    CollectGF=zeros(numel(q)/2,1);
	    ForceNum=9;
	    if ForceNum>0
	        CollectGF=zeros(numel(q)/2,ForceNum);
	    end
	    SubGF=zeros(numel(q)/2,1);
	    for FCount=1:ForceNum
	%SWITCHCASE_
			switch FCount
			    case 1
			        SubFrame=[0];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ForceFrameVecJac=ForceFVJac_1(t,q,p,u,s,SubSubs);
			        SubJac=ForceJac_1(t,q,p,u,s,SubSubs);
			        SubEff=ForceEff_1(t,q,p,u,s,SubSubs);
			        if SubFrame(1)~=0
			            JacSubs=Jac_Global(:,:,SubFrame);
			            for sfCount=1:numel(SubFrame)
			                SubJac=SubJac+ForceFrameVecJac(:,:,sfCount)*JacSubs(:,:,sfCount);
			            end
			        end
			        SubGF=SubJac.'*SubEff;
			    case 2
			        SubFrame=[0];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ForceFrameVecJac=ForceFVJac_2(t,q,p,u,s,SubSubs);
			        SubJac=ForceJac_2(t,q,p,u,s,SubSubs);
			        SubEff=ForceEff_2(t,q,p,u,s,SubSubs);
			        if SubFrame(1)~=0
			            JacSubs=Jac_Global(:,:,SubFrame);
			            for sfCount=1:numel(SubFrame)
			                SubJac=SubJac+ForceFrameVecJac(:,:,sfCount)*JacSubs(:,:,sfCount);
			            end
			        end
			        SubGF=SubJac.'*SubEff;
			    case 3
			        SubFrame=[7];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[7];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_3(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			    case 4
			        SubFrame=[21];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[21];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_4(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			    case 5
			        SubFrame=[22];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[22];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_5(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			    case 6
			        SubFrame=[23];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[23];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_6(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			    case 7
			        SubFrame=[24];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[24];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_7(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			    case 8
			        SubFrame=[25];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[25];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_8(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			    case 9
			        SubFrame=[26];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[26];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_9(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			end
	        CollectGF(:,FCount)=SubGF;
	    end
		function out1 = ForceJac_1(t,in2,in3,in4,in5,SUBSVECTOR__)
		%FORCEJAC_1
		%    OUT1 = FORCEJAC_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:50
		outSize=[3,9];
		elementRow=[1,2,3];
		elementCol=[1,2,3];
		elementList=[1.0,1.0,1.0];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = ForceJac_2(t,in2,in3,in4,in5,SUBSVECTOR__)
		%FORCEJAC_2
		%    OUT1 = FORCEJAC_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:51
		outSize=[6,9];
		elementRow=[1,2,3,4,5,6];
		elementCol=[4,5,6,7,8,9];
		elementList=[1.0,1.0,1.0,1.0,1.0,1.0];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = ForceJac_3(t,in2,in3,in4,in5,in6)
		%FORCEJAC_3
		%    OUT1 = FORCEJAC_3(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:52
		out1 = 0.0;
		end
		function out1 = ForceJac_4(t,in2,in3,in4,in5,in6)
		%FORCEJAC_4
		%    OUT1 = FORCEJAC_4(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:53
		out1 = 0.0;
		end
		function out1 = ForceJac_5(t,in2,in3,in4,in5,in6)
		%FORCEJAC_5
		%    OUT1 = FORCEJAC_5(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:54
		out1 = 0.0;
		end
		function out1 = ForceJac_6(t,in2,in3,in4,in5,in6)
		%FORCEJAC_6
		%    OUT1 = FORCEJAC_6(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:55
		out1 = 0.0;
		end
		function out1 = ForceJac_7(t,in2,in3,in4,in5,in6)
		%FORCEJAC_7
		%    OUT1 = FORCEJAC_7(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:55
		out1 = 0.0;
		end
		function out1 = ForceJac_8(t,in2,in3,in4,in5,in6)
		%FORCEJAC_8
		%    OUT1 = FORCEJAC_8(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:56
		out1 = 0.0;
		end
		function out1 = ForceJac_9(t,in2,in3,in4,in5,in6)
		%FORCEJAC_9
		%    OUT1 = FORCEJAC_9(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:57
		out1 = 0.0;
		end
		function out1 = ForceEff_1(t,in2,in3,in4,in5,SUBSVECTOR__)
		%FORCEEFF_1
		%    OUT1 = FORCEEFF_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:50
		b1act = in3(4,:);
		q1dev__dt_1_ = in2(10,:);
		q1flex__dt_1_ = in2(11,:);
		q1sup__dt_1_ = in2(12,:);
		out1 = [-b1act.*q1dev__dt_1_;-b1act.*q1flex__dt_1_;-b1act.*q1sup__dt_1_];
		end
		function out1 = ForceEff_2(t,in2,in3,in4,in5,SUBSVECTOR__)
		%FORCEEFF_2
		%    OUT1 = FORCEEFF_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:51
		b2act = in3(5,:);
		q2w1__dt_1_ = in2(13,:);
		q2w2__dt_1_ = in2(14,:);
		q2w3__dt_1_ = in2(15,:);
		q2w4__dt_1_ = in2(16,:);
		q2w5__dt_1_ = in2(17,:);
		q2w6__dt_1_ = in2(18,:);
		out1 = [-b2act.*q2w1__dt_1_;-b2act.*q2w2__dt_1_;-b2act.*q2w3__dt_1_;-b2act.*q2w4__dt_1_;-b2act.*q2w5__dt_1_;-b2act.*q2w6__dt_1_];
		end
		function out1 = ForceEff_3(t,in2,in3,in4,in5,in6)
		%FORCEEFF_3
		%    OUT1 = FORCEEFF_3(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:52
		gUncertain = in3(2,:);
		m13 = in3(63,:);
		out1 = [0.0;0.0;-gUncertain.*m13];
		end
		function out1 = ForceEff_4(t,in2,in3,in4,in5,in6)
		%FORCEEFF_4
		%    OUT1 = FORCEEFF_4(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:53
		g = in3(1,:);
		m22 = in3(70,:);
		out1 = [0.0;0.0;-g.*m22];
		end
		function out1 = ForceEff_5(t,in2,in3,in4,in5,in6)
		%FORCEEFF_5
		%    OUT1 = FORCEEFF_5(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:54
		g = in3(1,:);
		m23 = in3(77,:);
		out1 = [0.0;0.0;-g.*m23];
		end
		function out1 = ForceEff_6(t,in2,in3,in4,in5,in6)
		%FORCEEFF_6
		%    OUT1 = FORCEEFF_6(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:55
		g = in3(1,:);
		m24 = in3(84,:);
		out1 = [0.0;0.0;-g.*m24];
		end
		function out1 = ForceEff_7(t,in2,in3,in4,in5,in6)
		%FORCEEFF_7
		%    OUT1 = FORCEEFF_7(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:56
		g = in3(1,:);
		m25 = in3(91,:);
		out1 = [0.0;0.0;-g.*m25];
		end
		function out1 = ForceEff_8(t,in2,in3,in4,in5,in6)
		%FORCEEFF_8
		%    OUT1 = FORCEEFF_8(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:57
		g = in3(1,:);
		m26 = in3(98,:);
		out1 = [0.0;0.0;-g.*m26];
		end
		function out1 = ForceEff_9(t,in2,in3,in4,in5,in6)
		%FORCEEFF_9
		%    OUT1 = FORCEEFF_9(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:58
		g = in3(1,:);
		m27 = in3(105,:);
		out1 = [0.0;0.0;-g.*m27];
		end
		function out1 = ForceFVJac_1(t,in2,in3,in4,in5,SUBSVECTOR__)
		%FORCEFVJAC_1
		%    OUT1 = FORCEFVJAC_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:51
		out1 = 0.0;
		end
		function out1 = ForceFVJac_2(t,in2,in3,in4,in5,SUBSVECTOR__)
		%FORCEFVJAC_2
		%    OUT1 = FORCEFVJAC_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:52
		out1 = 0.0;
		end
		function out1 = ForceFVJac_3(t,in2,in3,in4,in5,in6)
		%FORCEFVJAC_3
		%    OUT1 = FORCEFVJAC_3(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:53
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = ForceFVJac_4(t,in2,in3,in4,in5,in6)
		%FORCEFVJAC_4
		%    OUT1 = FORCEFVJAC_4(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:53
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = ForceFVJac_5(t,in2,in3,in4,in5,in6)
		%FORCEFVJAC_5
		%    OUT1 = FORCEFVJAC_5(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:54
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = ForceFVJac_6(t,in2,in3,in4,in5,in6)
		%FORCEFVJAC_6
		%    OUT1 = FORCEFVJAC_6(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:55
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = ForceFVJac_7(t,in2,in3,in4,in5,in6)
		%FORCEFVJAC_7
		%    OUT1 = FORCEFVJAC_7(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:56
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = ForceFVJac_8(t,in2,in3,in4,in5,in6)
		%FORCEFVJAC_8
		%    OUT1 = FORCEFVJAC_8(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:57
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = ForceFVJac_9(t,in2,in3,in4,in5,in6)
		%FORCEFVJAC_9
		%    OUT1 = FORCEFVJAC_9(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:58
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
	end

	function [CollectGF,CollectInputJac]=Input_TAWE(t,q,p,u,s,TF_Global,TransDis_Global,Vel_Global,Jac_Global,Quat_Global)
	    CollectGF=zeros(numel(q)/2,1);
	    CollectInputJac=zeros(numel(q)/2,numel(u));
	    ForceNum=9;
	    if ForceNum>0
	        CollectGF=zeros(numel(q)/2,ForceNum);
	        CollectInputJac=zeros(numel(q)/2,numel(u),ForceNum);
	    end
	    SubGF=zeros(numel(q)/2,1);
	    SubInputJac=zeros(numel(q)/2,numel(u));
	    for FCount=1:ForceNum
	%SWITCHCASE_
			switch FCount
			    case 1
			        SubFrame=[8];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[8];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_1(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_1(t,q,p,u,s,SubSubs);
			    case 2
			        SubFrame=[9];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[9];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_2(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_2(t,q,p,u,s,SubSubs);
			    case 3
			        SubFrame=[10];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[10];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_3(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_3(t,q,p,u,s,SubSubs);
			    case 4
			        SubFrame=[11];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[11];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_4(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_4(t,q,p,u,s,SubSubs);
			    case 5
			        SubFrame=[7];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-7];
			        RefFrame=[3];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_5(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_5(t,q,p,u,s,SubSubs);
			    case 6
			        SubFrame=[7];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-7];
			        RefFrame=[4];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_6(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_6(t,q,p,u,s,SubSubs);
			    case 7
			        SubFrame=[21];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-21];
			        RefFrame=[14];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_7(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_7(t,q,p,u,s,SubSubs);
			    case 8
			        SubFrame=[21];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-21];
			        RefFrame=[15];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_8(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_8(t,q,p,u,s,SubSubs);
			    case 9
			        SubFrame=[22];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-22];
			        RefFrame=[15];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_9(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_9(t,q,p,u,s,SubSubs);
			end
	        CollectGF(:,FCount)=SubGF;
	        CollectInputJac(:,:,FCount)=SubInputJac;
	    end
		function out1 = ForceJac_1(t,in2,in3,in4,in5,in6)
		%FORCEJAC_1
		%    OUT1 = FORCEJAC_1(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:58
		out1 = 0.0;
		end
		function out1 = ForceJac_2(t,in2,in3,in4,in5,in6)
		%FORCEJAC_2
		%    OUT1 = FORCEJAC_2(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:59
		out1 = 0.0;
		end
		function out1 = ForceJac_3(t,in2,in3,in4,in5,in6)
		%FORCEJAC_3
		%    OUT1 = FORCEJAC_3(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:01
		out1 = 0.0;
		end
		function out1 = ForceJac_4(t,in2,in3,in4,in5,in6)
		%FORCEJAC_4
		%    OUT1 = FORCEJAC_4(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:02
		out1 = 0.0;
		end
		function out1 = ForceJac_5(t,in2,in3,in4,in5,in6)
		%FORCEJAC_5
		%    OUT1 = FORCEJAC_5(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:03
		out1 = 0.0;
		end
		function out1 = ForceJac_6(t,in2,in3,in4,in5,in6)
		%FORCEJAC_6
		%    OUT1 = FORCEJAC_6(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:05
		out1 = 0.0;
		end
		function out1 = ForceJac_7(t,in2,in3,in4,in5,in6)
		%FORCEJAC_7
		%    OUT1 = FORCEJAC_7(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:06
		out1 = 0.0;
		end
		function out1 = ForceJac_8(t,in2,in3,in4,in5,in6)
		%FORCEJAC_8
		%    OUT1 = FORCEJAC_8(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:08
		out1 = 0.0;
		end
		function out1 = ForceJac_9(t,in2,in3,in4,in5,in6)
		%FORCEJAC_9
		%    OUT1 = FORCEJAC_9(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:09
		out1 = 0.0;
		end
		function out1 = InputEff_1(t,in2,in3,in4,in5,in6)
		%INPUTEFF_1
		%    OUT1 = INPUTEFF_1(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:59
		g = in3(1,:);
		handmpt1 = in4(5,:);
		out1 = [0.0;0.0;-g.*handmpt1];
		end
		function out1 = InputEff_2(t,in2,in3,in4,in5,in6)
		%INPUTEFF_2
		%    OUT1 = INPUTEFF_2(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:00
		g = in3(1,:);
		handmpt2 = in4(6,:);
		out1 = [0.0;0.0;-g.*handmpt2];
		end
		function out1 = InputEff_3(t,in2,in3,in4,in5,in6)
		%INPUTEFF_3
		%    OUT1 = INPUTEFF_3(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:01
		g = in3(1,:);
		handmpt3 = in4(7,:);
		out1 = [0.0;0.0;-g.*handmpt3];
		end
		function out1 = InputEff_4(t,in2,in3,in4,in5,in6)
		%INPUTEFF_4
		%    OUT1 = INPUTEFF_4(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:02
		g = in3(1,:);
		handmpt4 = in4(8,:);
		out1 = [0.0;0.0;-g.*handmpt4];
		end
		function out1 = InputEff_5(t,in2,in3,in4,in5,in6)
		%INPUTEFF_5
		%    OUT1 = INPUTEFF_5(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:03
		u1dev = in4(1,:);
		out1 = [0.0;0.0;u1dev];
		end
		function out1 = InputEff_6(t,in2,in3,in4,in5,in6)
		%INPUTEFF_6
		%    OUT1 = INPUTEFF_6(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:05
		u1flex = in4(2,:);
		out1 = [u1flex;0.0;0.0];
		end
		function out1 = InputEff_7(t,in2,in3,in4,in5,in6)
		%INPUTEFF_7
		%    OUT1 = INPUTEFF_7(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:07
		u21 = in4(3,:);
		out1 = [0.0;0.0;u21];
		end
		function out1 = InputEff_8(t,in2,in3,in4,in5,in6)
		%INPUTEFF_8
		%    OUT1 = INPUTEFF_8(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:08
		u22 = in4(4,:);
		out1 = [-u22;0.0;0.0];
		end
		function out1 = InputEff_9(t,in2,in3,in4,in5,in6)
		%INPUTEFF_9
		%    OUT1 = INPUTEFF_9(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:10
		u22 = in4(4,:);
		out1 = [u22;0.0;0.0];
		end
		function out1 = InputEffJac_1(t,in2,in3,in4,in5,in6)
		%INPUTEFFJAC_1
		%    OUT1 = INPUTEFFJAC_1(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:59
		g = in3(1,:);
		outSize=[3,8];
		elementRow=[3];
		elementCol=[5];
		elementList=[-g];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_2(t,in2,in3,in4,in5,in6)
		%INPUTEFFJAC_2
		%    OUT1 = INPUTEFFJAC_2(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:00
		g = in3(1,:);
		outSize=[3,8];
		elementRow=[3];
		elementCol=[6];
		elementList=[-g];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_3(t,in2,in3,in4,in5,in6)
		%INPUTEFFJAC_3
		%    OUT1 = INPUTEFFJAC_3(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:01
		g = in3(1,:);
		outSize=[3,8];
		elementRow=[3];
		elementCol=[7];
		elementList=[-g];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_4(t,in2,in3,in4,in5,in6)
		%INPUTEFFJAC_4
		%    OUT1 = INPUTEFFJAC_4(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:03
		g = in3(1,:);
		outSize=[3,8];
		elementRow=[3];
		elementCol=[8];
		elementList=[-g];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_5(t,in2,in3,in4,in5,in6)
		%INPUTEFFJAC_5
		%    OUT1 = INPUTEFFJAC_5(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:04
		rf1IMU1x = in3(24,:);
		rf1IMU1y = in3(25,:);
		rf1IMU1z = in3(26,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		rf1Devz = in3(29,:);
		t2 = (cos(rf1Devx));
		t3 = (sin(rf1IMU1x));
		t4 = (sin(rf1IMU1z));
		t5 = (cos(rf1IMU1x));
		t6 = (cos(rf1IMU1z));
		t7 = (sin(rf1IMU1y));
		t8 = (sin(rf1Devx));
		t9 = (sin(rf1Devz));
		t10 = (cos(rf1Devz));
		t11 = (sin(rf1Devy));
		t12 = (t8.*t10);
		t26 = (t2.*t9.*t11);
		t13 = (t12-t26);
		t16 = (rf1IMU1x./2.0);
		t14 = (cos(t16));
		t18 = (rf1IMU1z./2.0);
		t15 = (cos(t18));
		t17 = (t14.^2);
		t19 = (t15.^2);
		t20 = (rf1IMU1y./2.0);
		t21 = (cos(rf1Devy));
		t22 = (cos(rf1IMU1y));
		t23 = (t8.*t9);
		t24 = (t2.*t10.*t11);
		t25 = (t23+t24);
		outSize=[3,8];
		elementRow=[1,2,3];
		elementCol=[1,1,1];
		elementList=[t13.*(t4.*t5-t3.*t6.*t7)+t2.*t21.*(t3.*t4+t5.*t6.*t7)+t6.*t22.*t25,-t13.*(t17.*-2.0-t19.*2.0+t17.*t19.*4.0+t14.*t15.*cos(t20).*sin(t16).*sin(t18).*sin(t20).*8.0+1.0)-t2.*t21.*(t3.*t6-t4.*t5.*t7)+t4.*t22.*t25,-t7.*t25-t3.*t13.*t22+t2.*t5.*t21.*t22];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_6(t,in2,in3,in4,in5,in6)
		%INPUTEFFJAC_6
		%    OUT1 = INPUTEFFJAC_6(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:06
		q1dev__dt_0_ = in2(1,:);
		q1sup__dt_0_ = in2(3,:);
		rf1IMU1x = in3(24,:);
		rf1IMU1y = in3(25,:);
		rf1IMU1z = in3(26,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		rf1Devz = in3(29,:);
		t4 = (rf1Devx./2.0);
		t2 = (cos(t4));
		t6 = (rf1Devz./2.0);
		t3 = (cos(t6));
		t5 = (t2.^2);
		t7 = (t3.^2);
		t8 = (rf1Devy./2.0);
		t9 = (cos(q1sup__dt_0_));
		t10 = (sin(rf1Devz));
		t11 = (sin(rf1IMU1x));
		t12 = (sin(rf1IMU1z));
		t13 = (cos(rf1IMU1x));
		t14 = (cos(rf1IMU1z));
		t15 = (sin(rf1IMU1y));
		t16 = (cos(q1dev__dt_0_));
		t17 = (sin(rf1Devy));
		t18 = (cos(rf1Devx));
		t19 = (cos(rf1Devy));
		t20 = (sin(q1sup__dt_0_));
		t21 = (sin(q1dev__dt_0_));
		t22 = (sin(rf1Devx));
		t23 = (cos(rf1Devz));
		t24 = (t22.*t23);
		t59 = (t10.*t17.*t18);
		t25 = (t24-t59);
		t26 = (t20.*t25);
		t27 = (t5.*t7.*4.0);
		t28 = (cos(t8));
		t29 = (sin(t4));
		t30 = (sin(t8));
		t31 = (sin(t6));
		t32 = (t2.*t3.*t28.*t29.*t30.*t31.*8.0);
		t60 = (t5.*2.0);
		t61 = (t7.*2.0);
		t33 = (t27+t32-t60-t61+1.0);
		t34 = (t9.*t21.*t33);
		t35 = (t9.*t10.*t16.*t19);
		t36 = (t26+t34+t35);
		t39 = (rf1IMU1x./2.0);
		t37 = (cos(t39));
		t41 = (rf1IMU1z./2.0);
		t38 = (cos(t41));
		t40 = (t37.^2);
		t42 = (t38.^2);
		t43 = (rf1IMU1y./2.0);
		t44 = (t9.*t16.*t17);
		t45 = (t18.*t19.*t20);
		t58 = (t9.*t19.*t21.*t22);
		t46 = (t44+t45-t58);
		t47 = (cos(rf1IMU1y));
		t48 = (t10.*t22);
		t49 = (t17.*t18.*t23);
		t50 = (t48+t49);
		t51 = (t20.*t50);
		t52 = (t10.*t18);
		t56 = (t17.*t22.*t23);
		t53 = (t52-t56);
		t54 = (t9.*t21.*t53);
		t57 = (t9.*t16.*t19.*t23);
		t55 = (t51+t54-t57);
		outSize=[3,8];
		elementRow=[1,2,3];
		elementCol=[2,2,2];
		elementList=[-t36.*(t12.*t13-t11.*t14.*t15)-t46.*(t11.*t12+t13.*t14.*t15)-t14.*t47.*t55,t46.*(t11.*t14-t12.*t13.*t15)+t36.*(t40.*-2.0-t42.*2.0+t40.*t42.*4.0+t37.*t38.*cos(t43).*sin(t39).*sin(t41).*sin(t43).*8.0+1.0)-t12.*t47.*t55,t15.*t55+t11.*t36.*t47-t13.*t46.*t47];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_7(t,in2,in3,in4,in5,in6)
		%INPUTEFFJAC_7
		%    OUT1 = INPUTEFFJAC_7(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:07
		rf1IMU1x = in3(24,:);
		rf1IMU1y = in3(25,:);
		rf1IMU1z = in3(26,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf2Medz = in3(35,:);
		t2 = (cos(rf2Medx));
		t3 = (sin(rf1IMU1x));
		t4 = (sin(rf1IMU1z));
		t5 = (cos(rf1IMU1x));
		t6 = (cos(rf1IMU1z));
		t7 = (sin(rf1IMU1y));
		t8 = (sin(rf2Medx));
		t9 = (sin(rf2Medz));
		t10 = (cos(rf2Medz));
		t11 = (sin(rf2Medy));
		t12 = (t8.*t10);
		t26 = (t2.*t9.*t11);
		t13 = (t12-t26);
		t16 = (rf1IMU1x./2.0);
		t14 = (cos(t16));
		t18 = (rf1IMU1z./2.0);
		t15 = (cos(t18));
		t17 = (t14.^2);
		t19 = (t15.^2);
		t20 = (rf1IMU1y./2.0);
		t21 = (cos(rf2Medy));
		t22 = (cos(rf1IMU1y));
		t23 = (t8.*t9);
		t24 = (t2.*t10.*t11);
		t25 = (t23+t24);
		outSize=[3,8];
		elementRow=[1,2,3];
		elementCol=[3,3,3];
		elementList=[t13.*(t4.*t5-t3.*t6.*t7)+t2.*t21.*(t3.*t4+t5.*t6.*t7)+t6.*t22.*t25,-t13.*(t17.*-2.0-t19.*2.0+t17.*t19.*4.0+t14.*t15.*cos(t20).*sin(t16).*sin(t18).*sin(t20).*8.0+1.0)-t2.*t21.*(t3.*t6-t4.*t5.*t7)+t4.*t22.*t25,-t7.*t25-t3.*t13.*t22+t2.*t5.*t21.*t22];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_8(t,in2,in3,in4,in5,in6)
		%INPUTEFFJAC_8
		%    OUT1 = INPUTEFFJAC_8(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:09
		q2w1__dt_0_ = in2(4,:);
		rf1IMU1x = in3(24,:);
		rf1IMU1y = in3(25,:);
		rf1IMU1z = in3(26,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf2Medz = in3(35,:);
		riw1 = in3(36,:);
		t4 = (rf2Medx./2.0);
		t2 = (cos(t4));
		t6 = (rf2Medz./2.0);
		t3 = (cos(t6));
		t5 = (t2.^2);
		t7 = (t3.^2);
		t8 = (rf2Medy./2.0);
		t9 = (q2w1__dt_0_+riw1);
		t10 = (cos(t9));
		t11 = (sin(t9));
		t12 = (cos(rf2Medy));
		t13 = (sin(rf1IMU1x));
		t14 = (sin(rf1IMU1z));
		t15 = (cos(rf1IMU1x));
		t16 = (cos(rf1IMU1z));
		t17 = (sin(rf1IMU1y));
		t18 = (sin(rf2Medz));
		t19 = (sin(rf2Medx));
		t20 = (sin(rf2Medy));
		t21 = (cos(rf2Medz));
		t22 = (t5.*t7.*4.0);
		t23 = (cos(t8));
		t24 = (sin(t4));
		t25 = (sin(t8));
		t26 = (sin(t6));
		t27 = (t2.*t3.*t23.*t24.*t25.*t26.*8.0);
		t50 = (t5.*2.0);
		t51 = (t7.*2.0);
		t28 = (t22+t27-t50-t51+1.0);
		t29 = (t11.*t28);
		t30 = (t10.*t12.*t18);
		t31 = (t29+t30);
		t34 = (rf1IMU1x./2.0);
		t32 = (cos(t34));
		t36 = (rf1IMU1z./2.0);
		t33 = (cos(t36));
		t35 = (t32.^2);
		t37 = (t33.^2);
		t38 = (rf1IMU1y./2.0);
		t39 = (t10.*t20);
		t49 = (t11.*t12.*t19);
		t40 = (t39-t49);
		t41 = (cos(rf1IMU1y));
		t42 = (cos(rf2Medx));
		t43 = (t18.*t42);
		t47 = (t19.*t20.*t21);
		t44 = (t43-t47);
		t45 = (t11.*t44);
		t48 = (t10.*t12.*t21);
		t46 = (t45-t48);
		outSize=[3,8];
		elementRow=[1,2,3];
		elementCol=[4,4,4];
		elementList=[t31.*(t14.*t15-t13.*t16.*t17)+t40.*(t13.*t14+t15.*t16.*t17)+t16.*t41.*t46,-t40.*(t13.*t16-t14.*t15.*t17)-t31.*(t35.*-2.0-t37.*2.0+t35.*t37.*4.0+t32.*t33.*cos(t38).*sin(t34).*sin(t36).*sin(t38).*8.0+1.0)+t14.*t41.*t46,-t17.*t46-t13.*t31.*t41+t15.*t40.*t41];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_9(t,in2,in3,in4,in5,in6)
		%INPUTEFFJAC_9
		%    OUT1 = INPUTEFFJAC_9(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:10
		q2w1__dt_0_ = in2(4,:);
		rf1IMU1x = in3(24,:);
		rf1IMU1y = in3(25,:);
		rf1IMU1z = in3(26,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf2Medz = in3(35,:);
		riw1 = in3(36,:);
		t4 = (rf2Medx./2.0);
		t2 = (cos(t4));
		t6 = (rf2Medz./2.0);
		t3 = (cos(t6));
		t5 = (t2.^2);
		t7 = (t3.^2);
		t8 = (rf2Medy./2.0);
		t9 = (q2w1__dt_0_+riw1);
		t10 = (cos(t9));
		t11 = (sin(t9));
		t12 = (cos(rf2Medy));
		t13 = (sin(rf1IMU1x));
		t14 = (sin(rf1IMU1z));
		t15 = (cos(rf1IMU1x));
		t16 = (cos(rf1IMU1z));
		t17 = (sin(rf1IMU1y));
		t18 = (sin(rf2Medz));
		t19 = (sin(rf2Medx));
		t20 = (sin(rf2Medy));
		t21 = (cos(rf2Medz));
		t22 = (t5.*t7.*4.0);
		t23 = (cos(t8));
		t24 = (sin(t4));
		t25 = (sin(t8));
		t26 = (sin(t6));
		t27 = (t2.*t3.*t23.*t24.*t25.*t26.*8.0);
		t50 = (t5.*2.0);
		t51 = (t7.*2.0);
		t28 = (t22+t27-t50-t51+1.0);
		t29 = (t11.*t28);
		t30 = (t10.*t12.*t18);
		t31 = (t29+t30);
		t34 = (rf1IMU1x./2.0);
		t32 = (cos(t34));
		t36 = (rf1IMU1z./2.0);
		t33 = (cos(t36));
		t35 = (t32.^2);
		t37 = (t33.^2);
		t38 = (rf1IMU1y./2.0);
		t39 = (t10.*t20);
		t49 = (t11.*t12.*t19);
		t40 = (t39-t49);
		t41 = (cos(rf1IMU1y));
		t42 = (cos(rf2Medx));
		t43 = (t18.*t42);
		t47 = (t19.*t20.*t21);
		t44 = (t43-t47);
		t45 = (t11.*t44);
		t48 = (t10.*t12.*t21);
		t46 = (t45-t48);
		outSize=[3,8];
		elementRow=[1,2,3];
		elementCol=[4,4,4];
		elementList=[-t31.*(t14.*t15-t13.*t16.*t17)-t40.*(t13.*t14+t15.*t16.*t17)-t16.*t41.*t46,t40.*(t13.*t16-t14.*t15.*t17)+t31.*(t35.*-2.0-t37.*2.0+t35.*t37.*4.0+t32.*t33.*cos(t38).*sin(t34).*sin(t36).*sin(t38).*8.0+1.0)-t14.*t41.*t46,t17.*t46+t13.*t31.*t41-t15.*t40.*t41];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputFVJac_1(t,in2,in3,in4,in5,in6)
		%INPUTFVJAC_1
		%    OUT1 = INPUTFVJAC_1(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:16:59
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_2(t,in2,in3,in4,in5,in6)
		%INPUTFVJAC_2
		%    OUT1 = INPUTFVJAC_2(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:00
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_3(t,in2,in3,in4,in5,in6)
		%INPUTFVJAC_3
		%    OUT1 = INPUTFVJAC_3(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:02
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_4(t,in2,in3,in4,in5,in6)
		%INPUTFVJAC_4
		%    OUT1 = INPUTFVJAC_4(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:03
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_5(t,in2,in3,in4,in5,in6)
		%INPUTFVJAC_5
		%    OUT1 = INPUTFVJAC_5(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:04
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
		function out1 = InputFVJac_6(t,in2,in3,in4,in5,in6)
		%INPUTFVJAC_6
		%    OUT1 = INPUTFVJAC_6(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:06
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
		function out1 = InputFVJac_7(t,in2,in3,in4,in5,in6)
		%INPUTFVJAC_7
		%    OUT1 = INPUTFVJAC_7(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:07
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
		function out1 = InputFVJac_8(t,in2,in3,in4,in5,in6)
		%INPUTFVJAC_8
		%    OUT1 = INPUTFVJAC_8(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:09
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
		function out1 = InputFVJac_9(t,in2,in3,in4,in5,in6)
		%INPUTFVJAC_9
		%    OUT1 = INPUTFVJAC_9(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:11
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
	end

	function [ConsJac,ConsCor,ConsGF]=Constraint_TAWE(t,q,p,u,s,TransDis_Global,Vel_Global,Cor_Global,Jac_Global,Quat_Global)
	%% Constraint dynamic property calculator (Toolbox Internal Use)
	%% 
	%% Robotics & Mechatronics Lab, Virginia Tech. (No Copyright Claimed)
	%% Author: Jiamin Wang; Revised: 28-Jan-2019
	    ConsNum=7;
	    if(ConsNum>0)
	        ConsJac=zeros(ConsNum,numel(q)/2);
	        ConsCor=zeros(ConsNum,1);
	        ConsGF=zeros(ConsNum,1);
	    else
	        ConsJac=zeros(1,numel(q)/2);
	        ConsCor=0;
	        ConsGF=0;
	    end
	    
	    SubConsJac=zeros(1,numel(q)/2);
	    SubConsCor=0;
	    SubConsGF=0;
	    for ConsCount=1:ConsNum
	        
	%SWITCHCASE_
			switch ConsCount
			    case 1
			        SubFrame=[0];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ConsFrameVecJac=ConsFVJac_1(t,q,p,u,s,SubSubs);
			        SubConsJac=ConsJac_1(t,q,p,u,s,SubSubs);
			        SubConsCor=ConsCor_1(t,q,p,u,s,SubSubs);
			        SubConsGF=ConsGF_1(t,q,p,u,s,SubSubs);
			        if SubFrame(1)~=0
			            JacSubs=Jac_Global(:,:,SubFrame);
			            CorSubs=Cor_Global(:,SubFrame);
			            for sfCount=1:numel(SubFrame)
			                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
			                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
			            end
			        end
			    case 2
			        SubFrame=[0];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ConsFrameVecJac=ConsFVJac_2(t,q,p,u,s,SubSubs);
			        SubConsJac=ConsJac_2(t,q,p,u,s,SubSubs);
			        SubConsCor=ConsCor_2(t,q,p,u,s,SubSubs);
			        SubConsGF=ConsGF_2(t,q,p,u,s,SubSubs);
			        if SubFrame(1)~=0
			            JacSubs=Jac_Global(:,:,SubFrame);
			            CorSubs=Cor_Global(:,SubFrame);
			            for sfCount=1:numel(SubFrame)
			                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
			                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
			            end
			        end
			    case 3
			        SubFrame=[0];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ConsFrameVecJac=ConsFVJac_3(t,q,p,u,s,SubSubs);
			        SubConsJac=ConsJac_3(t,q,p,u,s,SubSubs);
			        SubConsCor=ConsCor_3(t,q,p,u,s,SubSubs);
			        SubConsGF=ConsGF_3(t,q,p,u,s,SubSubs);
			        if SubFrame(1)~=0
			            JacSubs=Jac_Global(:,:,SubFrame);
			            CorSubs=Cor_Global(:,SubFrame);
			            for sfCount=1:numel(SubFrame)
			                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
			                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
			            end
			        end
			    case 4
			        SubFrame=[5  20];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ConsFrameVecJac=ConsFVJac_4(t,q,p,u,s,SubSubs);
			        SubConsJac=ConsJac_4(t,q,p,u,s,SubSubs);
			        SubConsCor=ConsCor_4(t,q,p,u,s,SubSubs);
			        SubConsGF=ConsGF_4(t,q,p,u,s,SubSubs);
			        if SubFrame(1)~=0
			            JacSubs=Jac_Global(:,:,SubFrame);
			            CorSubs=Cor_Global(:,SubFrame);
			            for sfCount=1:numel(SubFrame)
			                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
			                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
			            end
			        end
			    case 5
			        SubFrame=[5  20];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ConsFrameVecJac=ConsFVJac_5(t,q,p,u,s,SubSubs);
			        SubConsJac=ConsJac_5(t,q,p,u,s,SubSubs);
			        SubConsCor=ConsCor_5(t,q,p,u,s,SubSubs);
			        SubConsGF=ConsGF_5(t,q,p,u,s,SubSubs);
			        if SubFrame(1)~=0
			            JacSubs=Jac_Global(:,:,SubFrame);
			            CorSubs=Cor_Global(:,SubFrame);
			            for sfCount=1:numel(SubFrame)
			                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
			                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
			            end
			        end
			    case 6
			        SubFrame=[5  20];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ConsFrameVecJac=ConsFVJac_6(t,q,p,u,s,SubSubs);
			        SubConsJac=ConsJac_6(t,q,p,u,s,SubSubs);
			        SubConsCor=ConsCor_6(t,q,p,u,s,SubSubs);
			        SubConsGF=ConsGF_6(t,q,p,u,s,SubSubs);
			        if SubFrame(1)~=0
			            JacSubs=Jac_Global(:,:,SubFrame);
			            CorSubs=Cor_Global(:,SubFrame);
			            for sfCount=1:numel(SubFrame)
			                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
			                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
			            end
			        end
			    case 7
			        SubFrame=[0];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ConsFrameVecJac=ConsFVJac_7(t,q,p,u,s,SubSubs);
			        SubConsJac=ConsJac_7(t,q,p,u,s,SubSubs);
			        SubConsCor=ConsCor_7(t,q,p,u,s,SubSubs);
			        SubConsGF=ConsGF_7(t,q,p,u,s,SubSubs);
			        if SubFrame(1)~=0
			            JacSubs=Jac_Global(:,:,SubFrame);
			            CorSubs=Cor_Global(:,SubFrame);
			            for sfCount=1:numel(SubFrame)
			                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
			                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
			            end
			        end
			end
	        
	        ConsJac(ConsCount,:)=SubConsJac;
	        ConsCor(ConsCount,:)=SubConsCor;
	        ConsGF(ConsCount,:)=SubConsGF;
	    end
		function out1 = ConsJac_1(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSJAC_1
		%    OUT1 = CONSJAC_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:15
		l21 = in3(19,:);
		l22 = in3(20,:);
		l23 = in3(21,:);
		l24 = in3(22,:);
		l25 = in3(23,:);
		l12x = in3(10,:);
		l12y = in3(11,:);
		l13x = in3(13,:);
		l12z = in3(12,:);
		l13y = in3(14,:);
		l13z = in3(15,:);
		q1dev__dt_0_ = in2(1,:);
		q1flex__dt_0_ = in2(2,:);
		q1sup__dt_0_ = in2(3,:);
		q2w1__dt_0_ = in2(4,:);
		q2w2__dt_0_ = in2(5,:);
		q2w3__dt_0_ = in2(6,:);
		q2w4__dt_0_ = in2(7,:);
		q2w5__dt_0_ = in2(8,:);
		q2w6__dt_0_ = in2(9,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf2Medz = in3(35,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		rf1Devz = in3(29,:);
		riw1 = in3(36,:);
		riw2 = in3(37,:);
		riw3 = in3(38,:);
		riw4 = in3(39,:);
		riw5 = in3(40,:);
		riw6 = in3(41,:);
		t2 = (q1dev__dt_0_./2.0);
		t3 = (cos(t2));
		t6 = (q1flex__dt_0_./2.0);
		t4 = (cos(t6));
		t5 = (sin(t2));
		t7 = (q1sup__dt_0_./2.0);
		t8 = (cos(t7));
		t9 = (sin(t6));
		t10 = (sin(t7));
		t11 = (cos(q1dev__dt_0_));
		t12 = (cos(rf1Devz));
		t13 = (cos(q1flex__dt_0_));
		t14 = (sin(q1dev__dt_0_));
		t15 = (sin(q1flex__dt_0_));
		t16 = (sin(q1sup__dt_0_));
		t17 = (cos(q1sup__dt_0_));
		t18 = (sin(rf1Devx));
		t19 = (sin(rf1Devz));
		t20 = (cos(rf1Devx));
		t21 = (sin(rf1Devy));
		t22 = (t19.*t20);
		t33 = (t12.*t18.*t21);
		t23 = (t22-t33);
		t24 = (t11.*t13);
		t25 = (t14.*t15.*t16);
		t26 = (t24+t25);
		t27 = (t3.^2);
		t28 = (t4.^2);
		t29 = (cos(rf1Devy));
		t30 = (t14.*t15);
		t31 = (t11.*t13.*t16);
		t32 = (t30+t31);
		t34 = (t18.*t19);
		t35 = (t12.*t20.*t21);
		t36 = (t34+t35);
		t37 = (q2w5__dt_0_+riw5);
		t38 = (q2w6__dt_0_+riw6);
		t39 = (q2w4__dt_0_+riw4);
		t40 = (cos(t38));
		t41 = (sin(t38));
		t42 = (q2w3__dt_0_+riw3);
		t43 = (cos(t39));
		t44 = (cos(t37));
		t45 = (l25+l13z);
		t46 = (t44.*t45);
		t47 = (sin(t37));
		t48 = (l13x.*t40);
		t57 = (l13y.*t41);
		t49 = (l24+t48-t57);
		t58 = (t47.*t49);
		t50 = (t46-t58);
		t51 = (sin(t39));
		t52 = (l13y.*t40);
		t53 = (l13x.*t41);
		t54 = (l23+t52+t53);
		t55 = (q2w2__dt_0_+riw2);
		t56 = (cos(t42));
		t59 = (t43.*t54);
		t69 = (t50.*t51);
		t60 = (l22+t59-t69);
		t61 = (sin(t42));
		t62 = (t43.*t50);
		t63 = (t51.*t54);
		t64 = (t62+t63);
		t65 = (q2w1__dt_0_+riw1);
		t66 = (cos(rf2Medz));
		t67 = (cos(t65));
		t68 = (sin(t55));
		t70 = (t60.*t61);
		t71 = (t56.*t64);
		t72 = (t70+t71);
		t73 = (t68.*t72);
		t74 = (cos(t55));
		t75 = (t56.*t60);
		t82 = (t61.*t64);
		t76 = (l21+t75-t82);
		t83 = (t74.*t76);
		t77 = (t73-t83);
		t78 = (sin(t65));
		t79 = (t45.*t47);
		t80 = (t44.*t49);
		t81 = (t79+t80);
		t84 = (sin(rf2Medx));
		t85 = (sin(rf2Medz));
		t86 = (cos(rf2Medx));
		t87 = (sin(rf2Medy));
		t88 = (t85.*t86);
		t97 = (t66.*t84.*t87);
		t89 = (t88-t97);
		t90 = (cos(rf2Medy));
		t91 = (t72.*t74);
		t92 = (t68.*t76);
		t93 = (t91+t92);
		t94 = (t84.*t85);
		t95 = (t66.*t86.*t87);
		t96 = (t94+t95);
		t98 = (t75-t82);
		t99 = (t68.*t98);
		t100 = (t91+t99);
		t103 = (t59-t69);
		t101 = (t61.*t103);
		t102 = (t71+t101);
		t104 = (t82-t56.*t103);
		t105 = (t68.*t104);
		t106 = (t43.*t56.*t81);
		t112 = (t51.*t61.*t81);
		t107 = (t106-t112);
		t108 = (t51.*t56.*t81);
		t109 = (t43.*t61.*t81);
		t110 = (t108+t109);
		t111 = (t74.*t110);
		t113 = (t68.*t107);
		t114 = (t111+t113);
		t115 = (t48-t57);
		t116 = (t52+t53);
		t117 = (t51.*t115);
		t118 = (t43.*t47.*t116);
		t119 = (t117+t118);
		t120 = (t43.*t115);
		t123 = (t47.*t51.*t116);
		t121 = (t120-t123);
		t122 = (t56.*t119);
		t124 = (t61.*t121);
		t125 = (t122+t124);
		t126 = (t56.*t121);
		t128 = (t61.*t119);
		t127 = (t126-t128);
		t129 = (t68.*t125);
		t130 = (t129-t74.*t127);
		outSize=[1,9];
		elementRow=[1,1,1,1,1,1,1,1,1];
		elementCol=[1,2,3,4,5,6,7,8,9];
		elementList=[-t23.*(l12z.*t32+l12y.*(t3.*t5.*2.0-t3.*t5.*t28.*4.0-t4.*t5.^2.*t8.*t9.*t10.*4.0+t4.*t8.*t9.*t10.*t27.*4.0)+l12x.*t11.*t17)-t12.*t29.*(-l12z.*(t11.*t15-t13.*t14.*t16)+l12y.*t26+l12x.*t14.*t17),t23.*(l12z.*t26-l12y.*(t4.*t9.*2.0-t4.*t9.*t27.*4.0-t3.*t5.*t8.*t9.^2.*t10.*4.0+t3.*t5.*t8.*t10.*t28.*4.0))+t36.*(l12y.*t13.*t17-l12z.*t15.*t17)+t12.*t29.*(l12z.*(t13.*t14-t11.*t15.*t16)+l12y.*t32),-t23.*(l12y.*(t3.*t4.*t5.*t8.^2.*t9.*4.0-t3.*t4.*t5.*t9.*t10.^2.*4.0)-l12x.*t14.*t16+l12z.*t13.*t14.*t17)-t36.*(l12x.*t17+l12z.*t13.*t16+l12y.*t15.*t16)+t12.*t29.*(-l12x.*t11.*t16+l12z.*t11.*t13.*t17+l12y.*t11.*t15.*t17),t89.*(t67.*t81+t77.*t78)-t66.*t90.*(t67.*t77-t78.*t81),t77.*t96-t67.*t89.*t93-t66.*t78.*t90.*t93,t96.*(t73-t74.*t98)-t67.*t89.*t100-t66.*t78.*t90.*t100,t96.*(t74.*(t82-t56.*(t59-t69))+t68.*t102)+t67.*t89.*(t105-t74.*t102)+t66.*t78.*t90.*(t105-t74.*(t71+t61.*(t59-t69))),t89.*(t50.*t78+t67.*t114)-t96.*(t68.*t110-t74.*t107)-t66.*t90.*(t50.*t67-t78.*t114),-t89.*(t67.*t130+t44.*t78.*t116)-t96.*(t68.*t127+t74.*t125)-t66.*t90.*(t78.*t130-t44.*t67.*t116)];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = ConsJac_2(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSJAC_2
		%    OUT1 = CONSJAC_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:18:03
		l21 = in3(19,:);
		l22 = in3(20,:);
		l23 = in3(21,:);
		l24 = in3(22,:);
		l25 = in3(23,:);
		l12x = in3(10,:);
		l12y = in3(11,:);
		l13x = in3(13,:);
		l12z = in3(12,:);
		l13y = in3(14,:);
		l13z = in3(15,:);
		q1dev__dt_0_ = in2(1,:);
		q1flex__dt_0_ = in2(2,:);
		q1sup__dt_0_ = in2(3,:);
		q2w1__dt_0_ = in2(4,:);
		q2w2__dt_0_ = in2(5,:);
		q2w3__dt_0_ = in2(6,:);
		q2w4__dt_0_ = in2(7,:);
		q2w5__dt_0_ = in2(8,:);
		q2w6__dt_0_ = in2(9,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf2Medz = in3(35,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		rf1Devz = in3(29,:);
		riw1 = in3(36,:);
		riw2 = in3(37,:);
		riw3 = in3(38,:);
		riw4 = in3(39,:);
		riw5 = in3(40,:);
		riw6 = in3(41,:);
		t2 = (q1dev__dt_0_./2.0);
		t3 = (cos(t2));
		t6 = (q1flex__dt_0_./2.0);
		t4 = (cos(t6));
		t5 = (sin(t2));
		t7 = (q1sup__dt_0_./2.0);
		t8 = (cos(t7));
		t9 = (sin(t6));
		t10 = (sin(t7));
		t11 = (cos(q1dev__dt_0_));
		t14 = (rf1Devx./2.0);
		t12 = (cos(t14));
		t16 = (rf1Devz./2.0);
		t13 = (cos(t16));
		t15 = (t12.^2);
		t17 = (t13.^2);
		t18 = (rf1Devy./2.0);
		t19 = (cos(q1flex__dt_0_));
		t20 = (sin(q1dev__dt_0_));
		t21 = (sin(q1flex__dt_0_));
		t22 = (sin(q1sup__dt_0_));
		t23 = (cos(q1sup__dt_0_));
		t24 = (sin(rf1Devz));
		t25 = (t11.*t19);
		t26 = (t20.*t21.*t22);
		t27 = (t25+t26);
		t28 = (t3.^2);
		t29 = (t4.^2);
		t30 = (t15.*t17.*4.0);
		t31 = (cos(t18));
		t32 = (sin(t14));
		t33 = (sin(t18));
		t34 = (sin(t16));
		t35 = (t12.*t13.*t31.*t32.*t33.*t34.*8.0);
		t41 = (t15.*2.0);
		t42 = (t17.*2.0);
		t36 = (t30+t35-t41-t42+1.0);
		t37 = (cos(rf1Devy));
		t38 = (t20.*t21);
		t39 = (t11.*t19.*t22);
		t40 = (t38+t39);
		t43 = (cos(rf1Devz));
		t44 = (sin(rf1Devx));
		t45 = (t43.*t44);
		t46 = (cos(rf1Devx));
		t47 = (sin(rf1Devy));
		t48 = (t45-t24.*t46.*t47);
		t49 = (q2w5__dt_0_+riw5);
		t50 = (q2w6__dt_0_+riw6);
		t51 = (q2w4__dt_0_+riw4);
		t52 = (cos(t50));
		t53 = (sin(t50));
		t54 = (q2w3__dt_0_+riw3);
		t55 = (cos(t51));
		t56 = (cos(t49));
		t57 = (l25+l13z);
		t58 = (t56.*t57);
		t59 = (sin(t49));
		t60 = (l13x.*t52);
		t69 = (l13y.*t53);
		t61 = (l24+t60-t69);
		t70 = (t59.*t61);
		t62 = (t58-t70);
		t63 = (sin(t51));
		t64 = (l13y.*t52);
		t65 = (l13x.*t53);
		t66 = (l23+t64+t65);
		t67 = (q2w2__dt_0_+riw2);
		t68 = (cos(t54));
		t71 = (t55.*t66);
		t87 = (t62.*t63);
		t72 = (l22+t71-t87);
		t73 = (sin(t54));
		t74 = (t55.*t62);
		t75 = (t63.*t66);
		t76 = (t74+t75);
		t77 = (q2w1__dt_0_+riw1);
		t80 = (rf2Medx./2.0);
		t78 = (cos(t80));
		t82 = (rf2Medz./2.0);
		t79 = (cos(t82));
		t81 = (t78.^2);
		t83 = (t79.^2);
		t84 = (rf2Medy./2.0);
		t85 = (cos(t77));
		t86 = (sin(t67));
		t88 = (t72.*t73);
		t89 = (t68.*t76);
		t90 = (t88+t89);
		t91 = (t86.*t90);
		t92 = (cos(t67));
		t93 = (t68.*t72);
		t100 = (t73.*t76);
		t94 = (l21+t93-t100);
		t101 = (t92.*t94);
		t95 = (t91-t101);
		t96 = (sin(t77));
		t97 = (t57.*t59);
		t98 = (t56.*t61);
		t99 = (t97+t98);
		t102 = (sin(rf2Medz));
		t103 = (t81.*t83.*4.0);
		t104 = (cos(t84));
		t105 = (sin(t80));
		t106 = (sin(t84));
		t107 = (sin(t82));
		t108 = (t78.*t79.*t104.*t105.*t106.*t107.*8.0);
		t121 = (t81.*2.0);
		t122 = (t83.*2.0);
		t109 = (t103+t108-t121-t122+1.0);
		t110 = (cos(rf2Medy));
		t111 = (t90.*t92);
		t112 = (t86.*t94);
		t113 = (t111+t112);
		t114 = (cos(rf2Medz));
		t115 = (sin(rf2Medx));
		t116 = (t114.*t115);
		t117 = (cos(rf2Medx));
		t118 = (sin(rf2Medy));
		t125 = (t102.*t117.*t118);
		t119 = (t116-t125);
		t120 = (t93-t100);
		t123 = (t86.*t120);
		t124 = (t111+t123);
		t126 = (t71-t87);
		t130 = (t68.*t126);
		t127 = (t100-t130);
		t128 = (t73.*(t71-t87));
		t129 = (t89+t128);
		t131 = (t86.*t127);
		t132 = (t63.*t68.*t99);
		t133 = (t55.*t73.*t99);
		t134 = (t132+t133);
		t135 = (t55.*t68.*t99);
		t138 = (t63.*t73.*t99);
		t136 = (t135-t138);
		t137 = (t92.*t134);
		t139 = (t86.*t136);
		t140 = (t137+t139);
		t141 = (t60-t69);
		t142 = (t64+t65);
		t143 = (t63.*t141);
		t144 = (t55.*t59.*t142);
		t145 = (t143+t144);
		t146 = (t55.*t141);
		t149 = (t59.*t63.*t142);
		t147 = (t146-t149);
		t148 = (t68.*t145);
		t150 = (t73.*t147);
		t151 = (t148+t150);
		t152 = (t68.*t147);
		t154 = (t73.*t145);
		t153 = (t152-t154);
		t155 = (t86.*t151);
		t156 = (t155-t92.*t153);
		outSize=[1,9];
		elementRow=[1,1,1,1,1,1,1,1,1];
		elementCol=[1,2,3,4,5,6,7,8,9];
		elementList=[t36.*(l12z.*t40+l12y.*(t3.*t5.*2.0-t3.*t5.*t29.*4.0-t4.*t5.^2.*t8.*t9.*t10.*4.0+t4.*t8.*t9.*t10.*t28.*4.0)+l12x.*t11.*t23)-t24.*t37.*(-l12z.*(t11.*t21-t19.*t20.*t22)+l12y.*t27+l12x.*t20.*t23),-t36.*(l12z.*t27-l12y.*(t4.*t9.*2.0-t4.*t9.*t28.*4.0-t3.*t5.*t8.*t9.^2.*t10.*4.0+t3.*t5.*t8.*t10.*t29.*4.0))-t48.*(l12y.*t19.*t23-l12z.*t21.*t23)+t24.*t37.*(l12z.*(t19.*t20-t11.*t21.*t22)+l12y.*t40),t36.*(l12y.*(t3.*t4.*t5.*t8.^2.*t9.*4.0-t3.*t4.*t5.*t9.*t10.^2.*4.0)-l12x.*t20.*t22+l12z.*t19.*t20.*t23)+t48.*(l12x.*t23+l12z.*t19.*t22+l12y.*t21.*t22)+t24.*t37.*(-l12x.*t11.*t22+l12z.*t11.*t19.*t23+l12y.*t11.*t21.*t23),-t109.*(t85.*t99+t95.*t96)-t102.*t110.*(t85.*t95-t96.*t99),-t95.*t119+t85.*t109.*t113-t96.*t102.*t110.*t113,-t119.*(t91-t92.*t120)+t85.*t109.*t124-t96.*t102.*t110.*t124,-t119.*(t92.*t127+t86.*(t89+t73.*t126))-t85.*t109.*(t131-t92.*t129)+t96.*t102.*t110.*(t131-t92.*t129),-t109.*(t62.*t96+t85.*t140)+t119.*(t86.*t134-t92.*t136)-t102.*t110.*(t62.*t85-t96.*t140),t109.*(t85.*t156+t56.*t96.*t142)+t119.*(t86.*t153+t92.*t151)-t102.*t110.*(t96.*t156-t56.*t85.*t142)];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = ConsJac_3(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSJAC_3
		%    OUT1 = CONSJAC_3(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:18:59
		l21 = in3(19,:);
		l22 = in3(20,:);
		l23 = in3(21,:);
		l24 = in3(22,:);
		l25 = in3(23,:);
		l12x = in3(10,:);
		l12y = in3(11,:);
		l13x = in3(13,:);
		l12z = in3(12,:);
		l13y = in3(14,:);
		l13z = in3(15,:);
		q1dev__dt_0_ = in2(1,:);
		q1flex__dt_0_ = in2(2,:);
		q1sup__dt_0_ = in2(3,:);
		q2w1__dt_0_ = in2(4,:);
		q2w2__dt_0_ = in2(5,:);
		q2w3__dt_0_ = in2(6,:);
		q2w4__dt_0_ = in2(7,:);
		q2w5__dt_0_ = in2(8,:);
		q2w6__dt_0_ = in2(9,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		riw1 = in3(36,:);
		riw2 = in3(37,:);
		riw3 = in3(38,:);
		riw4 = in3(39,:);
		riw5 = in3(40,:);
		riw6 = in3(41,:);
		t2 = (cos(q1dev__dt_0_));
		t3 = (sin(q1flex__dt_0_));
		t4 = (cos(q1flex__dt_0_));
		t5 = (sin(q1dev__dt_0_));
		t6 = (sin(q1sup__dt_0_));
		t7 = (q1dev__dt_0_./2.0);
		t8 = (cos(t7));
		t11 = (q1flex__dt_0_./2.0);
		t9 = (cos(t11));
		t10 = (sin(t7));
		t12 = (q1sup__dt_0_./2.0);
		t13 = (cos(t12));
		t14 = (sin(t11));
		t15 = (sin(t12));
		t16 = (cos(q1sup__dt_0_));
		t17 = (sin(rf1Devy));
		t18 = (t3.*t5);
		t19 = (t2.*t4.*t6);
		t20 = (t18+t19);
		t21 = (cos(rf1Devy));
		t22 = (sin(rf1Devx));
		t23 = (t2.*t4);
		t24 = (t3.*t5.*t6);
		t25 = (t23+t24);
		t26 = (t8.^2);
		t27 = (t9.^2);
		t28 = (cos(rf1Devx));
		t29 = (q2w5__dt_0_+riw5);
		t30 = (q2w6__dt_0_+riw6);
		t31 = (q2w4__dt_0_+riw4);
		t32 = (cos(t30));
		t33 = (sin(t30));
		t34 = (q2w3__dt_0_+riw3);
		t35 = (cos(t31));
		t36 = (cos(t29));
		t37 = (l25+l13z);
		t38 = (t36.*t37);
		t39 = (sin(t29));
		t40 = (l13x.*t32);
		t49 = (l13y.*t33);
		t41 = (l24+t40-t49);
		t50 = (t39.*t41);
		t42 = (t38-t50);
		t43 = (sin(t31));
		t44 = (l13y.*t32);
		t45 = (l13x.*t33);
		t46 = (l23+t44+t45);
		t47 = (q2w2__dt_0_+riw2);
		t48 = (cos(t34));
		t51 = (t35.*t46);
		t60 = (t42.*t43);
		t52 = (l22+t51-t60);
		t53 = (sin(t34));
		t54 = (t35.*t42);
		t55 = (t43.*t46);
		t56 = (t54+t55);
		t57 = (q2w1__dt_0_+riw1);
		t58 = (sin(t57));
		t59 = (sin(t47));
		t61 = (t52.*t53);
		t62 = (t48.*t56);
		t63 = (t61+t62);
		t64 = (t59.*t63);
		t65 = (cos(t47));
		t66 = (t48.*t52);
		t74 = (t53.*t56);
		t67 = (l21+t66-t74);
		t76 = (t65.*t67);
		t68 = (t64-t76);
		t69 = (cos(t57));
		t70 = (t37.*t39);
		t71 = (t36.*t41);
		t72 = (t70+t71);
		t73 = (sin(rf2Medy));
		t75 = (cos(rf2Medy));
		t77 = (sin(rf2Medx));
		t78 = (t63.*t65);
		t79 = (t59.*t67);
		t80 = (t78+t79);
		t81 = (cos(rf2Medx));
		t82 = (t66-t74);
		t83 = (t59.*t82);
		t84 = (t78+t83);
		t85 = (t51-t60);
		t89 = (t48.*t85);
		t86 = (t74-t89);
		t87 = (t53.*(t51-t60));
		t88 = (t62+t87);
		t90 = (t59.*t86);
		t91 = (t43.*t48.*t72);
		t92 = (t35.*t53.*t72);
		t93 = (t91+t92);
		t94 = (t65.*t93);
		t95 = (t35.*t48.*t72);
		t99 = (t43.*t53.*t72);
		t96 = (t95-t99);
		t97 = (t59.*t96);
		t98 = (t94+t97);
		t100 = (t40-t49);
		t101 = (t44+t45);
		t102 = (t43.*t100);
		t103 = (t35.*t39.*t101);
		t104 = (t102+t103);
		t105 = (t35.*t100);
		t107 = (t39.*t43.*t101);
		t106 = (t105-t107);
		t108 = (t48.*t106);
		t115 = (t53.*t104);
		t109 = (t108-t115);
		t110 = (t65.*t109);
		t111 = (t48.*t104);
		t112 = (t53.*t106);
		t113 = (t111+t112);
		t114 = (t110-t59.*t113);
		outSize=[1,9];
		elementRow=[1,1,1,1,1,1,1,1,1];
		elementCol=[1,2,3,4,5,6,7,8,9];
		elementList=[t17.*(-l12z.*(t2.*t3-t4.*t5.*t6)+l12y.*t25+l12x.*t5.*t16)+t21.*t22.*(l12z.*t20+l12y.*(t8.*t10.*2.0-t8.*t10.*t27.*4.0-t9.*t10.^2.*t13.*t14.*t15.*4.0+t9.*t13.*t14.*t15.*t26.*4.0)+l12x.*t2.*t16),-t17.*(l12z.*(t4.*t5-t2.*t3.*t6)+l12y.*t20)+t21.*t28.*(l12y.*t4.*t16-l12z.*t3.*t16)-t21.*t22.*(l12z.*t25-l12y.*(t9.*t14.*2.0-t9.*t14.*t26.*4.0-t8.*t10.*t13.*t14.^2.*t15.*4.0+t8.*t10.*t13.*t15.*t27.*4.0)),-t17.*(-l12x.*t2.*t6+l12y.*t2.*t3.*t16+l12z.*t2.*t4.*t16)-t21.*t28.*(l12x.*t16+l12y.*t3.*t6+l12z.*t4.*t6)+t21.*t22.*(l12y.*(t8.*t9.*t10.*t13.^2.*t14.*4.0-t8.*t9.*t10.*t14.*t15.^2.*4.0)-l12x.*t5.*t6+l12z.*t4.*t5.*t16),-t73.*(t58.*t72-t68.*t69)-t75.*t77.*(t58.*t68+t69.*t72),t58.*t73.*t80+t68.*t75.*t81+t69.*t75.*t77.*t80,t75.*t81.*(t64-t65.*t82)+t58.*t73.*t84+t69.*t75.*t77.*t84,-t58.*t73.*(t90-t65.*(t62+t53.*t85))+t75.*t81.*(t59.*t88+t65.*t86)-t69.*t75.*t77.*(t90-t65.*t88),t73.*(t42.*t69-t58.*t98)-t75.*t77.*(t42.*t58+t69.*t98)-t75.*t81.*(t59.*t93-t65.*t96),-t73.*(t58.*t114+t36.*t69.*t101)-t75.*t81.*(t59.*t109+t65.*t113)-t75.*t77.*(t69.*t114-t36.*t58.*t101)];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = ConsJac_4(t,in2,in3,in4,in5,in6)
		%CONSJAC_4
		%    OUT1 = CONSJAC_4(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:37
		outSize=[1,9];
		elementRow=[];
		elementCol=[];
		elementList=[];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = ConsJac_5(t,in2,in3,in4,in5,in6)
		%CONSJAC_5
		%    OUT1 = CONSJAC_5(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:40
		outSize=[1,9];
		elementRow=[];
		elementCol=[];
		elementList=[];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = ConsJac_6(t,in2,in3,in4,in5,in6)
		%CONSJAC_6
		%    OUT1 = CONSJAC_6(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:43
		outSize=[1,9];
		elementRow=[];
		elementCol=[];
		elementList=[];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = ConsJac_7(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSJAC_7
		%    OUT1 = CONSJAC_7(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:46
		q1dev__dt_0_ = in2(1,:);
		q1flex__dt_0_ = in2(2,:);
		q1sup__dt_0_ = in2(3,:);
		xx3__dt_0_ = in5(1,:);
		zz3__dt_0_ = in5(3,:);
		t2 = (q1flex__dt_0_./2.0);
		t3 = (q1dev__dt_0_./2.0);
		t4 = (q1sup__dt_0_./2.0);
		t5 = (cos(q1sup__dt_0_));
		t6 = (cos(t4));
		t7 = (t6.^2);
		t8 = (sin(t3));
		t9 = (sin(t2));
		t10 = (sin(t4));
		t11 = (cos(t3));
		t12 = (cos(t2));
		t13 = (t11.^2);
		t14 = (t12.^2);
		t15 = (t7.^2);
		t16 = (xx3__dt_0_./2.0);
		t17 = (zz3__dt_0_./2.0);
		t18 = (t7.*t13.*2.0);
		t19 = (t7.*t14.*2.0);
		t20 = (t8.*t9.*t10.*2.0);
		t21 = (t6.*t7.*t11.*t12.*2.0);
		t34 = (t7.*3.0);
		t35 = (t7.*t8.*t9.*t10.*2.0);
		t22 = (-t13-t14+t15+t18+t19+t20+t21-t34-t35+2.0);
		t23 = (1.0./t22);
		t24 = (cos(t16));
		t25 = (sin(t17));
		t26 = (cos(q1dev__dt_0_));
		t27 = ((t5.*t26)./2.0);
		t28 = (cos(q1flex__dt_0_));
		t29 = ((t5.*t28)./2.0);
		t30 = (t15.*2.0);
		t31 = ((t8.*t9.*t10)./2.0);
		t32 = (t6.*t7.*t11.*t12.*3.0);
		t33 = (-t5+t27+t29+t30+t31+t32-t6.*t11.*t12-t5.*t8.*t9.*t10.*(3.0./2.0));
		t36 = (cos(t17));
		t37 = (sin(t16));
		t38 = (sin(q1sup__dt_0_));
		t39 = (t38./2.0);
		t40 = ((t10.*t11.*t12)./2.0);
		t41 = (t6.*t8.*t9);
		t42 = ((t5.*t10.*t11.*t12)./2.0);
		t43 = (t39+t40+t41+t42-t6.*t7.*t8.*t9);
		t44 = (sin(q1flex__dt_0_));
		t45 = (sin(q1dev__dt_0_));
		outSize=[1,9];
		elementRow=[1,1,1];
		elementCol=[1,2,3];
		elementList=[t6.*t9.*t11.*(-1.0./2.0)+(t8.*t10.*t12)./2.0-(t23.*t24.*t25.*t43)./4.0+(t23.*t33.*t36.*t37)./4.0,t6.*t8.*t12.*(-1.0./2.0)+(t9.*t10.*t11)./2.0+(t23.*t24.*t25.*t33)./4.0-(t23.*t36.*t37.*t43)./4.0,t31-(t6.*t11.*t12)./2.0-(t23.*t24.*t25.*(t45./2.0+(t38.*t44)./2.0+t6.*t8.*t12.*2.0+t9.*t10.*t11.*(3.0./2.0)-t6.*t7.*t8.*t12+(t5.*t9.*t10.*t11)./2.0))./4.0-(t23.*t36.*t37.*(t44./2.0+(t38.*t45)./2.0+t6.*t9.*t11.*2.0+t8.*t10.*t12.*(3.0./2.0)-t6.*t7.*t9.*t11+(t5.*t8.*t10.*t12)./2.0))./4.0];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = ConsCor_1(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSCOR_1
		%    OUT1 = CONSCOR_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:50
		l21 = in3(19,:);
		l22 = in3(20,:);
		l23 = in3(21,:);
		l24 = in3(22,:);
		l25 = in3(23,:);
		l12x = in3(10,:);
		l12y = in3(11,:);
		l13x = in3(13,:);
		l12z = in3(12,:);
		l13y = in3(14,:);
		l13z = in3(15,:);
		q1dev__dt_0_ = in2(1,:);
		q1dev__dt_1_ = in2(10,:);
		q1flex__dt_0_ = in2(2,:);
		q1flex__dt_1_ = in2(11,:);
		q1sup__dt_0_ = in2(3,:);
		q1sup__dt_1_ = in2(12,:);
		q2w1__dt_0_ = in2(4,:);
		q2w1__dt_1_ = in2(13,:);
		q2w2__dt_0_ = in2(5,:);
		q2w2__dt_1_ = in2(14,:);
		q2w3__dt_0_ = in2(6,:);
		q2w3__dt_1_ = in2(15,:);
		q2w4__dt_0_ = in2(7,:);
		q2w4__dt_1_ = in2(16,:);
		q2w5__dt_0_ = in2(8,:);
		q2w5__dt_1_ = in2(17,:);
		q2w6__dt_0_ = in2(9,:);
		q2w6__dt_1_ = in2(18,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf2Medz = in3(35,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		rf1Devz = in3(29,:);
		riw1 = in3(36,:);
		riw2 = in3(37,:);
		riw3 = in3(38,:);
		riw4 = in3(39,:);
		riw5 = in3(40,:);
		riw6 = in3(41,:);
		t2 = sin(q1dev__dt_0_);
		t3 = q1dev__dt_0_./2.0;
		t4 = q1flex__dt_0_./2.0;
		t5 = q1sup__dt_0_./2.0;
		t6 = sin(rf1Devx);
		t7 = sin(rf1Devz);
		t8 = cos(rf1Devx);
		t9 = cos(rf1Devz);
		t10 = sin(rf1Devy);
		t11 = sin(q1sup__dt_0_);
		t12 = cos(q1flex__dt_0_);
		t13 = cos(q1sup__dt_0_);
		t14 = cos(q1dev__dt_0_);
		t15 = sin(q1flex__dt_0_);
		t16 = t7.*t8;
		t31 = t6.*t9.*t10;
		t17 = t16-t31;
		t18 = cos(t4);
		t19 = cos(t5);
		t20 = sin(t3);
		t21 = sin(t4);
		t22 = cos(t3);
		t23 = sin(t5);
		t24 = t20.^2;
		t25 = t23.^2;
		t26 = t22.^2;
		t27 = t19.^2;
		t28 = cos(rf1Devy);
		t29 = t21.^2;
		t30 = t18.^2;
		t32 = t6.*t7;
		t33 = t8.*t9.*t10;
		t34 = t32+t33;
		t35 = q2w6__dt_0_+riw6;
		t36 = q2w4__dt_0_+riw4;
		t37 = q2w5__dt_0_+riw5;
		t38 = cos(t37);
		t39 = q2w3__dt_0_+riw3;
		t40 = cos(t35);
		t41 = l13y.*t40;
		t42 = sin(t35);
		t43 = l13x.*t42;
		t44 = t41+t43;
		t45 = q2w2__dt_0_+riw2;
		t46 = cos(t39);
		t47 = cos(t36);
		t48 = sin(t39);
		t49 = sin(t36);
		t50 = q2w1__dt_0_+riw1;
		t51 = cos(t45);
		t52 = t38.*t44.*t46.*t47;
		t65 = t38.*t44.*t48.*t49;
		t53 = t52-t65;
		t54 = sin(t45);
		t55 = t38.*t44.*t46.*t49;
		t56 = t38.*t44.*t47.*t48;
		t57 = t55+t56;
		t58 = sin(rf2Medx);
		t59 = sin(rf2Medz);
		t60 = cos(rf2Medx);
		t61 = cos(rf2Medz);
		t62 = sin(rf2Medy);
		t63 = sin(t50);
		t64 = t51.*t57;
		t66 = t53.*t54;
		t67 = t64+t66;
		t68 = cos(t50);
		t69 = sin(t37);
		t70 = l13x.*t40;
		t73 = l13y.*t42;
		t71 = t70-t73;
		t72 = t44.*t49;
		t83 = t47.*t69.*t71;
		t74 = t72-t83;
		t75 = t44.*t47;
		t76 = t49.*t69.*t71;
		t77 = t75+t76;
		t78 = t59.*t60;
		t99 = t58.*t61.*t62;
		t79 = t78-t99;
		t80 = t58.*t59;
		t81 = t60.*t61.*t62;
		t82 = t80+t81;
		t84 = t46.*t74;
		t85 = t48.*t77;
		t86 = t84+t85;
		t87 = t48.*t74;
		t90 = t46.*t77;
		t88 = t87-t90;
		t89 = cos(rf2Medy);
		t91 = t51.*t88;
		t92 = t54.*t86;
		t93 = t91+t92;
		t94 = t49.*t71;
		t95 = t44.*t47.*t69;
		t96 = t94+t95;
		t97 = t47.*t71;
		t100 = t44.*t49.*t69;
		t98 = t97-t100;
		t101 = t46.*t98;
		t108 = t48.*t96;
		t102 = t101-t108;
		t103 = t51.*t102;
		t104 = t46.*t96;
		t105 = t48.*t98;
		t106 = t104+t105;
		t109 = t54.*t106;
		t107 = t103-t109;
		t110 = t51.*t106;
		t111 = t54.*t102;
		t112 = t110+t111;
		t113 = t82.*t107;
		t114 = t68.*t79.*t112;
		t115 = t61.*t63.*t89.*t112;
		t116 = t113+t114+t115;
		t117 = l25+l13z;
		t118 = t69.*t117;
		t119 = l24+t70-t73;
		t120 = t38.*t119;
		t121 = t118+t120;
		t122 = t46.*t47.*t121;
		t127 = t48.*t49.*t121;
		t123 = t122-t127;
		t124 = t46.*t49.*t121;
		t125 = t47.*t48.*t121;
		t126 = t124+t125;
		t128 = t51.*t123;
		t166 = t54.*t126;
		t129 = t128-t166;
		t130 = t38.*t117;
		t133 = t69.*t119;
		t131 = t130-t133;
		t132 = l23+t41+t43;
		t134 = t47.*t132;
		t139 = t49.*t131;
		t135 = l22+t134-t139;
		t136 = t47.*t131;
		t137 = t49.*t132;
		t138 = t136+t137;
		t140 = t46.*t135;
		t145 = t48.*t138;
		t141 = t140-t145;
		t142 = t48.*t135;
		t143 = t46.*t138;
		t144 = t142+t143;
		t146 = t51.*t141;
		t151 = t54.*t144;
		t147 = t146-t151;
		t148 = t51.*t144;
		t149 = t54.*t141;
		t150 = t148+t149;
		t152 = t68.*t79.*t147;
		t153 = t61.*t63.*t89.*t147;
		t171 = t82.*t150;
		t154 = t152+t153-t171;
		t155 = t134-t139;
		t157 = t46.*t155;
		t156 = t145-t157;
		t158 = t51.*t156;
		t159 = t48.*(t134-t139);
		t160 = t143+t159;
		t161 = t54.*t160;
		t162 = t158+t161;
		t163 = t51.*t126;
		t164 = t54.*t123;
		t165 = t163+t164;
		t167 = t68.*t79.*t129;
		t168 = t61.*t63.*t89.*t129;
		t211 = t82.*t165;
		t169 = t167+t168-t211;
		t170 = q2w5__dt_1_.*t169.*2.0;
		t172 = l21+t140-t145;
		t173 = t54.*t172;
		t174 = t148+t173;
		t210 = t51.*t172;
		t175 = t151-t210;
		t176 = t54.*t156;
		t177 = t68.*t79.*t162;
		t178 = t61.*t63.*t89.*t162;
		t179 = t14.*t15;
		t246 = t2.*t11.*t12;
		t180 = t179-t246;
		t181 = l12z.*t12.*t13;
		t182 = l12y.*t13.*t15;
		t183 = t2.*t12;
		t250 = t11.*t14.*t15;
		t184 = t183-t250;
		t185 = t20.*t22.*t25.*t29.*2.0;
		t186 = t20.*t22.*t27.*t30.*2.0;
		t187 = t185+t186-t20.*t22.*t25.*t30.*2.0-t20.*t22.*t27.*t29.*2.0;
		t188 = l12y.*t187;
		t189 = t188-l12z.*t2.*t13.*t15;
		t190 = t17.*t189;
		t191 = l12y.*t11.*t12;
		t192 = t191-l12z.*t11.*t15;
		t193 = t34.*t192;
		t194 = l12y.*t12.*t13.*t14;
		t195 = t194-l12z.*t13.*t14.*t15;
		t196 = t190+t193-t9.*t28.*t195;
		t197 = t63.*t79.*t174;
		t198 = t197-t61.*t68.*t89.*t174;
		t199 = t63.*t79.*t150;
		t200 = t199-t61.*t68.*t89.*t150;
		t201 = t63.*t107;
		t202 = t38.*t44.*t68;
		t203 = t201+t202;
		t204 = t79.*t203;
		t205 = t68.*t107;
		t206 = t205-t38.*t44.*t63;
		t207 = t204-t61.*t89.*t206;
		t209 = t51.*t160;
		t208 = t176-t209;
		t213 = t82.*t208;
		t212 = t177+t178-t213;
		t214 = q2w4__dt_1_.*t212.*2.0;
		t215 = t63.*t79.*t208;
		t216 = t215-t61.*t68.*t89.*t208;
		t217 = t67.*t68;
		t218 = t217-t44.*t63.*t69;
		t219 = t79.*t218;
		t220 = t51.*t53;
		t221 = t220-t54.*t57;
		t222 = t82.*t221;
		t223 = t63.*t67;
		t224 = t44.*t68.*t69;
		t225 = t223+t224;
		t226 = t61.*t89.*t225;
		t227 = t219+t222+t226;
		t228 = t63.*t121;
		t229 = t46.*t47.*t131;
		t235 = t48.*t49.*t131;
		t230 = t229-t235;
		t231 = t46.*t49.*t131;
		t232 = t47.*t48.*t131;
		t233 = t231+t232;
		t234 = t51.*t233;
		t236 = t54.*t230;
		t237 = t234+t236;
		t238 = t68.*t121;
		t239 = t68.*t131;
		t240 = t239-t63.*t165;
		t241 = t79.*t240;
		t242 = t63.*t131;
		t243 = t68.*t165;
		t244 = t242+t243;
		t245 = t61.*t89.*t244;
		t247 = t26.*t30.*2.0;
		t248 = t18.*t19.*t20.*t21.*t22.*t23.*8.0;
		t249 = l12x.*t2.*t13;
		t251 = t2.*t15;
		t252 = t11.*t12.*t14;
		t253 = t251+t252;
		t254 = l12x.*t13.*t14;
		t255 = l12y.*t184;
		t256 = l12z.*t184;
		t257 = t19.*t23.*t24.*t29.*2.0;
		t258 = t18.*t20.*t21.*t22.*4.0;
		t259 = t19.*t23.*t26.*t30.*2.0;
		t260 = t257+t258+t259-t19.*t23.*t24.*t30.*2.0-t19.*t23.*t26.*t29.*2.0;
		t261 = l12y.*t260;
		t262 = t256+t261;
		t263 = t17.*t262;
		t264 = l12y.*t180;
		t265 = t12.*t14;
		t266 = t2.*t11.*t15;
		t267 = t265+t266;
		t268 = l12z.*t267;
		t269 = t264+t268;
		t270 = t263-t9.*t28.*t269;
		t271 = t18.*t21.*t24.*t25.*2.0;
		t272 = t18.*t21.*t26.*t27.*2.0;
		t273 = t271+t272-t18.*t21.*t24.*t27.*2.0-t18.*t21.*t25.*t26.*2.0;
		t274 = l12y.*t273;
		t275 = l12z.*t12.*t13.*t14;
		t276 = t274+t275-l12x.*t11.*t14;
		t277 = t17.*t276;
		t278 = l12z.*t2.*t12.*t13;
		t279 = l12y.*t2.*t13.*t15;
		t280 = t278+t279-l12x.*t2.*t11;
		t281 = t9.*t28.*t280;
		t282 = t277+t281;
		t283 = l12z.*t180;
		out1 = q2w4__dt_1_.*(t170+t214-q2w6__dt_1_.*t116.*2.0-q2w1__dt_1_.*t216.*2.0+q2w2__dt_1_.*t212.*2.0+q2w3__dt_1_.*t212.*2.0)-q1dev__dt_1_.*(q1flex__dt_1_.*t270.*2.0+q1sup__dt_1_.*t282.*2.0-q1dev__dt_1_.*(t17.*(t249-t283+l12y.*(t24-t26+t247+t248-t24.*t30.*2.0))-t9.*t28.*(t254-t255+l12z.*t253)).*2.0)+q2w2__dt_1_.*(t170+t214-q2w6__dt_1_.*t116.*2.0-q2w3__dt_1_.*t154.*2.0+q2w1__dt_1_.*t198.*2.0+q2w2__dt_1_.*(t82.*t174+t68.*t79.*t175+t61.*t63.*t89.*t175).*2.0)-q2w6__dt_1_.*(q2w2__dt_1_.*t116.*2.0+q2w3__dt_1_.*t116.*2.0+q2w4__dt_1_.*t116.*2.0+q2w1__dt_1_.*t207.*2.0+q2w5__dt_1_.*t227.*2.0-q2w6__dt_1_.*(t79.*(t68.*t93-t38.*t63.*t71)+t82.*(t51.*t86-t54.*t88)+t61.*t89.*(t63.*t93+t38.*t68.*t71)).*2.0)+q2w3__dt_1_.*(t170-q2w2__dt_1_.*t154.*2.0-q2w6__dt_1_.*t116.*2.0-q2w3__dt_1_.*t154.*2.0+q2w1__dt_1_.*t200.*2.0+q2w4__dt_1_.*(t177+t178-t82.*(t176-t51.*(t143+t48.*t155))).*2.0)+q2w5__dt_1_.*(q2w1__dt_1_.*(t241+t245).*2.0+q2w5__dt_1_.*(t82.*(t51.*t230-t54.*t233)-t79.*(t228-t68.*t237)+t61.*t89.*(t238+t63.*t237)).*2.0+q2w2__dt_1_.*t169.*2.0+q2w3__dt_1_.*t169.*2.0+q2w4__dt_1_.*t169.*2.0-q2w6__dt_1_.*t227.*2.0)-q1flex__dt_1_.*(q1dev__dt_1_.*t270.*2.0+q1sup__dt_1_.*t196.*2.0+q1flex__dt_1_.*(t34.*(t181+t182)+t17.*(t283-l12y.*(t29-t30+t247+t248-t26.*t29.*2.0))-t9.*t28.*(t255-l12z.*t253)).*2.0)+q2w1__dt_1_.*(q2w5__dt_1_.*(t241+t245).*2.0+q2w2__dt_1_.*t198.*2.0+q2w3__dt_1_.*t200.*2.0-q2w4__dt_1_.*t216.*2.0-q2w6__dt_1_.*t207.*2.0-q2w1__dt_1_.*(t79.*(t228-t68.*t175)-t61.*t89.*(t238+t63.*t175)).*2.0)-q1sup__dt_1_.*(q1dev__dt_1_.*t282.*2.0+q1flex__dt_1_.*t196.*2.0+q1sup__dt_1_.*(-t17.*(t249+l12z.*t2.*t11.*t12+l12y.*t18.*t19.*t20.*t21.*t22.*t23.*8.0)+t34.*(t181+t182-l12x.*t11)+t9.*t28.*(t254+l12z.*t11.*t12.*t14+l12y.*t11.*t14.*t15)).*2.0);
		end
		function out1 = ConsCor_2(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSCOR_2
		%    OUT1 = CONSCOR_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:18:46
		l21 = in3(19,:);
		l22 = in3(20,:);
		l23 = in3(21,:);
		l24 = in3(22,:);
		l25 = in3(23,:);
		l12x = in3(10,:);
		l12y = in3(11,:);
		l13x = in3(13,:);
		l12z = in3(12,:);
		l13y = in3(14,:);
		l13z = in3(15,:);
		q1dev__dt_0_ = in2(1,:);
		q1dev__dt_1_ = in2(10,:);
		q1flex__dt_0_ = in2(2,:);
		q1flex__dt_1_ = in2(11,:);
		q1sup__dt_0_ = in2(3,:);
		q1sup__dt_1_ = in2(12,:);
		q2w1__dt_0_ = in2(4,:);
		q2w1__dt_1_ = in2(13,:);
		q2w2__dt_0_ = in2(5,:);
		q2w2__dt_1_ = in2(14,:);
		q2w3__dt_0_ = in2(6,:);
		q2w3__dt_1_ = in2(15,:);
		q2w4__dt_0_ = in2(7,:);
		q2w4__dt_1_ = in2(16,:);
		q2w5__dt_0_ = in2(8,:);
		q2w5__dt_1_ = in2(17,:);
		q2w6__dt_0_ = in2(9,:);
		q2w6__dt_1_ = in2(18,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf2Medz = in3(35,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		rf1Devz = in3(29,:);
		riw1 = in3(36,:);
		riw2 = in3(37,:);
		riw3 = in3(38,:);
		riw4 = in3(39,:);
		riw5 = in3(40,:);
		riw6 = in3(41,:);
		t8 = q1dev__dt_0_./2.0;
		t2 = cos(t8);
		t5 = q1flex__dt_0_./2.0;
		t3 = sin(t5);
		t4 = q1sup__dt_0_./2.0;
		t6 = cos(t5);
		t7 = cos(t4);
		t9 = sin(t8);
		t10 = sin(t4);
		t11 = t9.^2;
		t12 = t3.^2;
		t13 = t2.^2;
		t14 = t6.^2;
		t17 = rf1Devx./2.0;
		t15 = cos(t17);
		t19 = rf1Devz./2.0;
		t16 = cos(t19);
		t18 = t15.^2;
		t20 = t16.^2;
		t21 = rf1Devy./2.0;
		t22 = cos(q1dev__dt_0_);
		t23 = sin(q1flex__dt_0_);
		t24 = cos(q1flex__dt_0_);
		t25 = sin(q1dev__dt_0_);
		t26 = sin(q1sup__dt_0_);
		t27 = cos(q1sup__dt_0_);
		t28 = sin(rf1Devz);
		t29 = t22.*t23;
		t226 = t24.*t25.*t26;
		t30 = t29-t226;
		t31 = t18.*t20.*4.0;
		t32 = cos(t21);
		t33 = sin(t17);
		t34 = sin(t21);
		t35 = sin(t19);
		t36 = t15.*t16.*t32.*t33.*t34.*t35.*8.0;
		t43 = t18.*2.0;
		t44 = t20.*2.0;
		t37 = t31+t36-t43-t44+1.0;
		t38 = cos(rf1Devy);
		t39 = t24.*t25;
		t229 = t22.*t23.*t26;
		t40 = t39-t229;
		t41 = t10.^2;
		t42 = t7.^2;
		t45 = cos(rf1Devz);
		t46 = sin(rf1Devx);
		t47 = t45.*t46;
		t48 = cos(rf1Devx);
		t49 = sin(rf1Devy);
		t283 = t28.*t48.*t49;
		t50 = t47-t283;
		t51 = q2w5__dt_0_+riw5;
		t52 = q2w6__dt_0_+riw6;
		t53 = q2w4__dt_0_+riw4;
		t54 = q2w3__dt_0_+riw3;
		t55 = sin(t51);
		t56 = l25+l13z;
		t57 = t55.*t56;
		t58 = cos(t51);
		t59 = cos(t52);
		t60 = l13x.*t59;
		t61 = sin(t52);
		t68 = l13y.*t61;
		t62 = l24+t60-t68;
		t63 = t58.*t62;
		t64 = t57+t63;
		t65 = q2w2__dt_0_+riw2;
		t66 = cos(t54);
		t67 = cos(t53);
		t69 = sin(t54);
		t70 = sin(t53);
		t71 = cos(t65);
		t72 = t64.*t66.*t67;
		t87 = t64.*t69.*t70;
		t73 = t72-t87;
		t74 = sin(t65);
		t75 = t64.*t66.*t70;
		t76 = t64.*t67.*t69;
		t77 = t75+t76;
		t80 = rf2Medx./2.0;
		t78 = cos(t80);
		t82 = rf2Medz./2.0;
		t79 = cos(t82);
		t81 = t78.^2;
		t83 = t79.^2;
		t84 = rf2Medy./2.0;
		t85 = q2w1__dt_0_+riw1;
		t86 = sin(rf2Medz);
		t88 = t71.*t73;
		t101 = t74.*t77;
		t89 = t88-t101;
		t90 = cos(rf2Medz);
		t91 = sin(rf2Medx);
		t92 = t90.*t91;
		t93 = cos(rf2Medx);
		t94 = sin(rf2Medy);
		t113 = t86.*t93.*t94;
		t95 = t92-t113;
		t96 = t71.*t77;
		t97 = t73.*t74;
		t98 = t96+t97;
		t99 = t95.*t98;
		t100 = cos(t85);
		t102 = t81.*t83.*4.0;
		t103 = cos(t84);
		t104 = sin(t80);
		t105 = sin(t84);
		t106 = sin(t82);
		t107 = t78.*t79.*t103.*t104.*t105.*t106.*8.0;
		t114 = t81.*2.0;
		t115 = t83.*2.0;
		t108 = t102+t107-t114-t115+1.0;
		t109 = sin(t85);
		t110 = cos(rf2Medy);
		t111 = t86.*t89.*t109.*t110;
		t116 = t89.*t100.*t108;
		t112 = t99+t111-t116;
		t117 = l13y.*t59;
		t118 = l13x.*t61;
		t119 = t117+t118;
		t120 = t58.*t66.*t67.*t119;
		t126 = t58.*t69.*t70.*t119;
		t121 = t120-t126;
		t122 = t58.*t66.*t70.*t119;
		t123 = t58.*t67.*t69.*t119;
		t124 = t122+t123;
		t125 = t71.*t124;
		t127 = t74.*t121;
		t128 = t125+t127;
		t129 = t56.*t58;
		t131 = t55.*t62;
		t130 = t129-t131;
		t132 = t66.*t70.*t130;
		t133 = t67.*t69.*t130;
		t134 = t132+t133;
		t135 = t66.*t67.*t130;
		t138 = t69.*t70.*t130;
		t136 = t135-t138;
		t137 = t71.*t134;
		t139 = t74.*t136;
		t140 = t137+t139;
		t141 = t100.*t128;
		t142 = t141-t55.*t109.*t119;
		t143 = t108.*t142;
		t144 = t71.*t121;
		t145 = t144-t74.*t124;
		t146 = t95.*t145;
		t147 = t109.*t128;
		t148 = t55.*t100.*t119;
		t149 = t147+t148;
		t150 = t143+t146-t86.*t110.*t149;
		t151 = t60-t68;
		t152 = t70.*t119;
		t157 = t55.*t67.*t151;
		t153 = t152-t157;
		t154 = t67.*t119;
		t155 = t55.*t70.*t151;
		t156 = t154+t155;
		t158 = t66.*t153;
		t159 = t69.*t156;
		t160 = t158+t159;
		t161 = t66.*t156;
		t163 = t69.*t153;
		t162 = t161-t163;
		t164 = t71.*t162;
		t165 = t164-t74.*t160;
		t166 = t70.*t151;
		t167 = t55.*t67.*t119;
		t168 = t166+t167;
		t169 = t67.*t151;
		t171 = t55.*t70.*t119;
		t170 = t169-t171;
		t172 = t66.*t170;
		t179 = t69.*t168;
		t173 = t172-t179;
		t174 = t71.*t173;
		t175 = t66.*t168;
		t176 = t69.*t170;
		t177 = t175+t176;
		t180 = t74.*t177;
		t178 = t174-t180;
		t181 = t71.*t177;
		t182 = t74.*t173;
		t183 = t181+t182;
		t184 = t95.*t178;
		t185 = t100.*t108.*t183;
		t187 = t86.*t109.*t110.*t183;
		t186 = t184+t185-t187;
		t188 = l23+t117+t118;
		t189 = t67.*t188;
		t194 = t70.*t130;
		t190 = l22+t189-t194;
		t191 = t67.*t130;
		t192 = t70.*t188;
		t193 = t191+t192;
		t195 = t66.*t190;
		t200 = t69.*t193;
		t196 = t195-t200;
		t197 = t69.*t190;
		t198 = t66.*t193;
		t199 = t197+t198;
		t201 = t71.*t196;
		t207 = t74.*t199;
		t202 = t201-t207;
		t203 = t71.*t199;
		t204 = t74.*t196;
		t205 = t203+t204;
		t206 = t95.*t205;
		t208 = t86.*t109.*t110.*t202;
		t218 = t100.*t108.*t202;
		t209 = t206+t208-t218;
		t210 = t189-t194;
		t212 = t66.*t210;
		t211 = t200-t212;
		t213 = t71.*t211;
		t214 = t69.*(t189-t194);
		t215 = t198+t214;
		t216 = t74.*t215;
		t217 = t213+t216;
		t219 = q2w3__dt_1_.*t209.*2.0;
		t220 = l21+t195-t200;
		t221 = t74.*t220;
		t222 = t203+t221;
		t265 = t71.*t220;
		t223 = t207-t265;
		t224 = t74.*t211;
		t225 = t86.*t109.*t110.*t217;
		t227 = t13.*t14.*2.0;
		t228 = t2.*t3.*t6.*t7.*t9.*t10.*8.0;
		t230 = t23.*t25;
		t231 = t22.*t24.*t26;
		t232 = t230+t231;
		t233 = l12z.*t40;
		t234 = t7.*t10.*t11.*t12.*2.0;
		t235 = t2.*t3.*t6.*t9.*4.0;
		t236 = t7.*t10.*t13.*t14.*2.0;
		t237 = t234+t235+t236-t7.*t10.*t11.*t14.*2.0-t7.*t10.*t12.*t13.*2.0;
		t238 = l12y.*t237;
		t239 = t233+t238;
		t240 = t37.*t239;
		t241 = l12y.*t30;
		t242 = t22.*t24;
		t243 = t23.*t25.*t26;
		t244 = t242+t243;
		t245 = l12z.*t244;
		t246 = t241+t245;
		t247 = t28.*t38.*t246;
		t248 = t240+t247;
		t249 = t108.*t109.*t222;
		t250 = t86.*t100.*t110.*t222;
		t251 = t249+t250;
		t252 = t108.*t109.*t205;
		t253 = t86.*t100.*t110.*t205;
		t254 = t252+t253;
		t255 = t109.*t178;
		t256 = t58.*t100.*t119;
		t257 = t255+t256;
		t258 = t108.*t257;
		t259 = t100.*t178;
		t260 = t259-t58.*t109.*t119;
		t261 = t86.*t110.*t260;
		t262 = t258+t261;
		t264 = t71.*t215;
		t263 = t224-t264;
		t266 = t64.*t109;
		t267 = t64.*t100;
		t268 = t98.*t109;
		t269 = t268-t100.*t130;
		t270 = t108.*t269;
		t271 = t109.*t130;
		t272 = t98.*t100;
		t273 = t271+t272;
		t274 = t86.*t110.*t273;
		t275 = t270+t274;
		t276 = t95.*t263;
		t278 = t100.*t108.*t217;
		t277 = t225+t276-t278;
		t279 = t108.*t109.*t263;
		t280 = t86.*t100.*t110.*t263;
		t281 = t279+t280;
		t282 = l12x.*t25.*t27;
		t284 = l12z.*t24.*t27;
		t285 = l12y.*t23.*t27;
		t286 = l12x.*t22.*t27;
		t287 = t3.*t6.*t11.*t41.*2.0;
		t288 = t3.*t6.*t13.*t42.*2.0;
		t289 = t287+t288-t3.*t6.*t11.*t42.*2.0-t3.*t6.*t13.*t41.*2.0;
		t290 = l12y.*t289;
		t291 = l12z.*t22.*t24.*t27;
		t292 = t290+t291-l12x.*t22.*t26;
		t293 = t37.*t292;
		t294 = l12z.*t24.*t25.*t27;
		t295 = l12y.*t23.*t25.*t27;
		t296 = t294+t295-l12x.*t25.*t26;
		t297 = t293-t28.*t38.*t296;
		t298 = t2.*t9.*t12.*t41.*2.0;
		t299 = t2.*t9.*t14.*t42.*2.0;
		t300 = t298+t299-t2.*t9.*t12.*t42.*2.0-t2.*t9.*t14.*t41.*2.0;
		t301 = l12y.*t300;
		t302 = t37.*(t301-l12z.*t23.*t25.*t27);
		t303 = l12y.*t24.*t26;
		t304 = t303-l12z.*t23.*t26;
		t305 = t50.*t304;
		t306 = l12y.*t22.*t24.*t27;
		t307 = t306-l12z.*t22.*t23.*t27;
		t308 = t28.*t38.*t307;
		t309 = t302+t305+t308;
		out1 = q2w1__dt_1_.*(q2w2__dt_1_.*t251.*-2.0-q2w3__dt_1_.*t254.*2.0+q2w4__dt_1_.*t281.*2.0+q2w6__dt_1_.*t262.*2.0+q2w5__dt_1_.*t275.*2.0+q2w1__dt_1_.*(t108.*(t266-t100.*t223)+t86.*t110.*(t267+t109.*t223)).*2.0)+q1dev__dt_1_.*(q1flex__dt_1_.*t248.*2.0+q1sup__dt_1_.*t297.*2.0-q1dev__dt_1_.*(t37.*(t282-l12z.*t30+l12y.*(t11-t13+t227+t228-t11.*t14.*2.0))+t28.*t38.*(t286-l12y.*t40+l12z.*t232)).*2.0)-q2w2__dt_1_.*(t219-q2w5__dt_1_.*t112.*2.0-q2w6__dt_1_.*t186.*2.0+q2w1__dt_1_.*t251.*2.0-q2w4__dt_1_.*t277.*2.0+q2w2__dt_1_.*(t95.*t222+t100.*t108.*t223-t86.*t109.*t110.*t223).*2.0)+q1flex__dt_1_.*(q1dev__dt_1_.*t248.*2.0+q1sup__dt_1_.*t309.*2.0+q1flex__dt_1_.*(t37.*(l12z.*t30-l12y.*(t12-t14+t227+t228-t12.*t13.*2.0))+t50.*(t284+t285)+t28.*t38.*(l12y.*t40-l12z.*t232)).*2.0)+q2w6__dt_1_.*(q2w5__dt_1_.*t150.*2.0+q2w2__dt_1_.*t186.*2.0+q2w3__dt_1_.*t186.*2.0+q2w4__dt_1_.*t186.*2.0+q2w1__dt_1_.*t262.*2.0-q2w6__dt_1_.*(-t108.*(t100.*t165+t58.*t109.*t151)+t95.*(t71.*t160+t74.*t162)+t86.*t110.*(t109.*t165-t58.*t100.*t151)).*2.0)-q2w3__dt_1_.*(t219-q2w5__dt_1_.*t112.*2.0+q2w2__dt_1_.*t209.*2.0-q2w6__dt_1_.*t186.*2.0+q2w1__dt_1_.*t254.*2.0-q2w4__dt_1_.*(t225+t95.*(t224-t71.*(t198+t69.*t210))-t100.*t108.*t217).*2.0)+q2w5__dt_1_.*(q2w5__dt_1_.*(-t95.*(t71.*t136-t74.*t134)+t108.*(t266-t100.*t140)+t86.*t110.*(t267+t109.*t140)).*2.0+q2w2__dt_1_.*t112.*2.0+q2w3__dt_1_.*t112.*2.0+q2w4__dt_1_.*t112.*2.0+q2w6__dt_1_.*t150.*2.0+q2w1__dt_1_.*t275.*2.0)+q2w4__dt_1_.*(q2w5__dt_1_.*t112.*2.0+q2w6__dt_1_.*t186.*2.0+q2w1__dt_1_.*t281.*2.0+q2w2__dt_1_.*t277.*2.0+q2w3__dt_1_.*t277.*2.0+q2w4__dt_1_.*t277.*2.0)+q1sup__dt_1_.*(q1dev__dt_1_.*t297.*2.0+q1flex__dt_1_.*t309.*2.0-q1sup__dt_1_.*(t37.*(t282+l12z.*t24.*t25.*t26+l12y.*t2.*t3.*t6.*t7.*t9.*t10.*8.0)-t50.*(t284+t285-l12x.*t26)+t28.*t38.*(t286+l12y.*t22.*t23.*t26+l12z.*t22.*t24.*t26)).*2.0);
		end
		function out1 = ConsCor_3(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSCOR_3
		%    OUT1 = CONSCOR_3(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:30
		l21 = in3(19,:);
		l22 = in3(20,:);
		l23 = in3(21,:);
		l24 = in3(22,:);
		l25 = in3(23,:);
		l12x = in3(10,:);
		l12y = in3(11,:);
		l13x = in3(13,:);
		l12z = in3(12,:);
		l13y = in3(14,:);
		l13z = in3(15,:);
		q1dev__dt_0_ = in2(1,:);
		q1dev__dt_1_ = in2(10,:);
		q1flex__dt_0_ = in2(2,:);
		q1flex__dt_1_ = in2(11,:);
		q1sup__dt_0_ = in2(3,:);
		q1sup__dt_1_ = in2(12,:);
		q2w1__dt_0_ = in2(4,:);
		q2w1__dt_1_ = in2(13,:);
		q2w2__dt_0_ = in2(5,:);
		q2w2__dt_1_ = in2(14,:);
		q2w3__dt_0_ = in2(6,:);
		q2w3__dt_1_ = in2(15,:);
		q2w4__dt_0_ = in2(7,:);
		q2w4__dt_1_ = in2(16,:);
		q2w5__dt_0_ = in2(8,:);
		q2w5__dt_1_ = in2(17,:);
		q2w6__dt_0_ = in2(9,:);
		q2w6__dt_1_ = in2(18,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		riw1 = in3(36,:);
		riw2 = in3(37,:);
		riw3 = in3(38,:);
		riw4 = in3(39,:);
		riw5 = in3(40,:);
		riw6 = in3(41,:);
		t2 = q2w6__dt_0_+riw6;
		t3 = q2w4__dt_0_+riw4;
		t4 = cos(t2);
		t5 = sin(t2);
		t6 = q2w3__dt_0_+riw3;
		t7 = sin(t3);
		t8 = l13x.*t4;
		t18 = l13y.*t5;
		t9 = t8-t18;
		t10 = cos(t3);
		t11 = q2w5__dt_0_+riw5;
		t12 = sin(t11);
		t13 = l13y.*t4;
		t14 = l13x.*t5;
		t15 = t13+t14;
		t16 = q2w2__dt_0_+riw2;
		t17 = cos(t6);
		t19 = t7.*t9;
		t20 = t10.*t12.*t15;
		t21 = t19+t20;
		t22 = sin(t6);
		t23 = t9.*t10;
		t28 = t7.*t12.*t15;
		t24 = t23-t28;
		t25 = q2w1__dt_0_+riw1;
		t26 = sin(t25);
		t27 = cos(t16);
		t29 = t17.*t24;
		t39 = t21.*t22;
		t30 = t29-t39;
		t31 = sin(t16);
		t32 = t17.*t21;
		t33 = t22.*t24;
		t34 = t32+t33;
		t41 = t27.*t30;
		t42 = t31.*t34;
		t35 = t41-t42;
		t36 = cos(t25);
		t37 = cos(t11);
		t38 = sin(rf2Medy);
		t40 = cos(rf2Medy);
		t43 = sin(rf2Medx);
		t44 = t27.*t34;
		t45 = t30.*t31;
		t46 = t44+t45;
		t47 = t26.*t38.*t46;
		t48 = cos(rf2Medx);
		t49 = t36.*t40.*t43.*t46;
		t51 = t35.*t40.*t48;
		t50 = t47+t49-t51;
		t52 = t10.*t15.*t17.*t37;
		t58 = t7.*t15.*t22.*t37;
		t53 = t52-t58;
		t54 = t7.*t15.*t17.*t37;
		t55 = t10.*t15.*t22.*t37;
		t56 = t54+t55;
		t57 = t27.*t56;
		t59 = t31.*t53;
		t60 = t57+t59;
		t61 = t7.*t15;
		t67 = t9.*t10.*t12;
		t62 = t61-t67;
		t63 = t10.*t15;
		t64 = t7.*t9.*t12;
		t65 = t63+t64;
		t66 = t17.*t65;
		t74 = t22.*t62;
		t68 = t66-t74;
		t69 = t27.*t68;
		t70 = t17.*t62;
		t71 = t22.*t65;
		t72 = t70+t71;
		t73 = t69-t31.*t72;
		t75 = l25+l13z;
		t76 = t37.*t75;
		t77 = l24+t8-t18;
		t80 = t12.*t77;
		t78 = t76-t80;
		t79 = l23+t13+t14;
		t81 = t10.*t78;
		t82 = t7.*t79;
		t83 = t81+t82;
		t84 = t7.*t78;
		t87 = t10.*t79;
		t85 = t84-t87;
		t86 = t17.*t83;
		t93 = t22.*t85;
		t88 = t86-t93;
		t89 = t17.*t85;
		t90 = t22.*t83;
		t91 = t89+t90;
		t92 = t27.*t91;
		t94 = t31.*t88;
		t95 = t92+t94;
		t96 = t12.*t75;
		t97 = t37.*t77;
		t98 = t96+t97;
		t99 = t7.*t17.*t98;
		t100 = t10.*t22.*t98;
		t101 = t99+t100;
		t102 = t10.*t17.*t98;
		t104 = t7.*t22.*t98;
		t103 = t102-t104;
		t105 = t27.*t103;
		t131 = t31.*t101;
		t106 = t105-t131;
		t107 = l22-t84+t87;
		t108 = t22.*t107;
		t109 = t86+t108;
		t111 = t17.*t107;
		t110 = t90-t111;
		t112 = t27.*t110;
		t113 = t31.*t109;
		t114 = t112+t113;
		t115 = t26.*t38.*t114;
		t116 = t27.*t109;
		t120 = t31.*t110;
		t117 = t116-t120;
		t118 = t36.*t40.*t43.*t114;
		t140 = t40.*t48.*t117;
		t119 = t115+t118-t140;
		t121 = l21-t90+t111;
		t179 = t27.*t121;
		t122 = t113-t179;
		t123 = t26.*t38.*t95;
		t124 = t27.*t88;
		t178 = t31.*t91;
		t125 = t124-t178;
		t126 = t36.*t40.*t43.*t95;
		t186 = t40.*t48.*t125;
		t127 = t123+t126-t186;
		t128 = q2w4__dt_1_.*t127.*2.0;
		t129 = t31.*t121;
		t130 = t116+t129;
		t132 = t26.*t38.*t106;
		t133 = t27.*t101;
		t134 = t31.*t103;
		t135 = t133+t134;
		t136 = t40.*t48.*t135;
		t137 = t36.*t40.*t43.*t106;
		t138 = t132+t136+t137;
		t139 = q2w5__dt_1_.*t138.*2.0;
		t141 = q2w3__dt_1_.*t119.*2.0;
		t142 = sin(q1dev__dt_0_);
		t143 = sin(q1flex__dt_0_);
		t144 = cos(q1dev__dt_0_);
		t145 = cos(q1flex__dt_0_);
		t146 = sin(q1sup__dt_0_);
		t151 = q1dev__dt_0_./2.0;
		t147 = cos(t151);
		t149 = q1flex__dt_0_./2.0;
		t148 = cos(t149);
		t150 = t148.^2;
		t152 = t147.^2;
		t153 = sin(t149);
		t154 = t153.^2;
		t155 = q1sup__dt_0_./2.0;
		t156 = cos(rf1Devy);
		t157 = cos(q1sup__dt_0_);
		t158 = sin(rf1Devy);
		t159 = sin(rf1Devx);
		t160 = sin(t151);
		t161 = sin(t155);
		t162 = cos(t155);
		t163 = t161.^2;
		t164 = t162.^2;
		t165 = cos(rf1Devx);
		t166 = t143.*t144;
		t191 = t142.*t145.*t146;
		t167 = t166-t191;
		t168 = t142.*t145;
		t187 = t143.*t144.*t146;
		t169 = t168-t187;
		t170 = t160.^2;
		t171 = t35.*t36;
		t172 = t171-t15.*t26.*t37;
		t173 = t38.*t172;
		t174 = t26.*t35;
		t175 = t15.*t36.*t37;
		t176 = t174+t175;
		t177 = t173-t40.*t43.*t176;
		t180 = t36.*t38.*t130;
		t181 = t180-t26.*t40.*t43.*t130;
		t182 = t36.*t38.*t117;
		t183 = t182-t26.*t40.*t43.*t117;
		t184 = t36.*t38.*t125;
		t185 = t184-t26.*t40.*t43.*t125;
		t188 = t142.*t143;
		t189 = t144.*t145.*t146;
		t190 = t188+t189;
		t192 = t150.*t152.*2.0;
		t193 = t147.*t148.*t153.*t160.*t161.*t162.*8.0;
		t194 = l12z.*t167;
		t195 = l12y.*t167;
		t196 = t144.*t145;
		t197 = t142.*t143.*t146;
		t198 = t196+t197;
		t199 = l12z.*t198;
		t200 = t195+t199;
		t201 = t158.*t200;
		t202 = l12z.*t169;
		t203 = t154.*t161.*t162.*t170.*2.0;
		t204 = t147.*t148.*t153.*t160.*4.0;
		t205 = t150.*t152.*t161.*t162.*2.0;
		t206 = t203+t204+t205-t152.*t154.*t161.*t162.*2.0-t150.*t161.*t162.*t170.*2.0;
		t207 = l12y.*t206;
		t208 = t202+t207;
		t209 = t201-t156.*t159.*t208;
		t210 = t36.*t98;
		t211 = t26.*t98;
		t212 = t7.*t17.*t78;
		t213 = t10.*t22.*t78;
		t214 = t212+t213;
		t215 = t27.*t214;
		t216 = t10.*t17.*t78;
		t220 = t7.*t22.*t78;
		t217 = t216-t220;
		t218 = t31.*t217;
		t219 = t215+t218;
		t221 = t26.*t78;
		t222 = t36.*t135;
		t223 = t221+t222;
		t224 = t38.*t223;
		t225 = t36.*t78;
		t226 = t40.*t43.*(t225-t26.*t135);
		t227 = t224+t226;
		t228 = t26.*t60;
		t229 = t12.*t15.*t36;
		t230 = t228+t229;
		t231 = t38.*t230;
		t232 = t27.*t53;
		t233 = t232-t31.*t56;
		t234 = t36.*t60;
		t235 = t234-t12.*t15.*t26;
		t236 = t40.*t43.*t235;
		t237 = t231+t236-t40.*t48.*t233;
		t238 = l12z.*t143.*t144.*t157;
		t239 = t238-l12y.*t144.*t145.*t157;
		t240 = t158.*t239;
		t241 = t147.*t150.*t160.*t163.*2.0;
		t242 = t147.*t154.*t160.*t164.*2.0;
		t243 = t241+t242-t147.*t150.*t160.*t164.*2.0-t147.*t154.*t160.*t163.*2.0;
		t244 = l12y.*t243;
		t245 = l12z.*t142.*t143.*t157;
		t246 = t244+t245;
		t247 = l12z.*t143.*t146;
		t248 = t247-l12y.*t145.*t146;
		t249 = t156.*t165.*t248;
		t250 = t240+t249-t156.*t159.*t246;
		t251 = l12x.*t144.*t157;
		t252 = l12z.*t145.*t157;
		t253 = l12y.*t143.*t157;
		t254 = l12x.*t142.*t157;
		t255 = l12z.*t142.*t145.*t157;
		t256 = l12y.*t142.*t143.*t157;
		t257 = t255+t256-l12x.*t142.*t146;
		t258 = t158.*t257;
		t259 = t148.*t153.*t164.*t170.*2.0;
		t260 = t148.*t152.*t153.*t163.*2.0;
		t261 = t259+t260-t148.*t152.*t153.*t164.*2.0-t148.*t153.*t163.*t170.*2.0;
		t262 = l12y.*t261;
		t263 = l12x.*t144.*t146;
		t264 = t262+t263-l12z.*t144.*t145.*t157;
		t265 = t258-t156.*t159.*t264;
		out1 = q2w1__dt_1_.*(q2w2__dt_1_.*t181.*2.0+q2w3__dt_1_.*t183.*2.0+q2w4__dt_1_.*t185.*2.0-q2w6__dt_1_.*t177.*2.0-q2w5__dt_1_.*t227.*2.0-q2w1__dt_1_.*(t38.*(t210+t26.*t122)-t40.*t43.*(t211-t36.*t122)).*2.0)-q2w4__dt_1_.*(t128+t139-q2w6__dt_1_.*t50.*2.0+q2w2__dt_1_.*t127.*2.0+q2w3__dt_1_.*t127.*2.0-q2w1__dt_1_.*t185.*2.0)-q2w5__dt_1_.*(q2w2__dt_1_.*t138.*2.0+q2w3__dt_1_.*t138.*2.0+q2w4__dt_1_.*t138.*2.0+q2w1__dt_1_.*t227.*2.0-q2w6__dt_1_.*t237.*2.0-q2w5__dt_1_.*(-t38.*(t210+t26.*t219)+t40.*t48.*(t27.*t217-t31.*t214)+t40.*t43.*(t211-t36.*t219)).*2.0)-q2w3__dt_1_.*(t128+t139+t141-q2w6__dt_1_.*t50.*2.0+q2w2__dt_1_.*t119.*2.0-q2w1__dt_1_.*t183.*2.0)+q1sup__dt_1_.*(q1dev__dt_1_.*t265.*2.0+q1flex__dt_1_.*t250.*2.0-q1sup__dt_1_.*(-t158.*(t251+l12y.*t143.*t144.*t146+l12z.*t144.*t145.*t146)+t156.*t159.*(t254+l12z.*t142.*t145.*t146+l12y.*t147.*t148.*t153.*t160.*t161.*t162.*8.0)+t156.*t165.*(t252+t253-l12x.*t146)).*2.0)+q1dev__dt_1_.*(q1flex__dt_1_.*t209.*-2.0+q1sup__dt_1_.*t265.*2.0+q1dev__dt_1_.*(t158.*(t251-l12y.*t169+l12z.*t190)-t156.*t159.*(-t194+t254+l12y.*(-t152+t170+t192+t193-t150.*t170.*2.0))).*2.0)+q2w6__dt_1_.*(q2w2__dt_1_.*t50.*2.0+q2w3__dt_1_.*t50.*2.0+q2w4__dt_1_.*t50.*2.0-q2w1__dt_1_.*t177.*2.0+q2w5__dt_1_.*t237.*2.0+q2w6__dt_1_.*(t38.*(t26.*t73-t9.*t36.*t37)+t40.*t48.*(t27.*t72+t31.*t68)+t40.*t43.*(t36.*t73+t9.*t26.*t37)).*2.0)-q1flex__dt_1_.*(q1flex__dt_1_.*(t158.*(l12y.*t169-l12z.*t190)+t156.*t165.*(t252+t253)-t156.*t159.*(t194-l12y.*(-t150+t154+t192+t193-t152.*t154.*2.0))).*2.0+q1dev__dt_1_.*t209.*2.0-q1sup__dt_1_.*t250.*2.0)-q2w2__dt_1_.*(t128+t139+t141-q2w6__dt_1_.*t50.*2.0-q2w1__dt_1_.*t181.*2.0+q2w2__dt_1_.*(t26.*t38.*t122-t40.*t48.*t130+t36.*t40.*t43.*t122).*2.0);
		end
		function out1 = ConsCor_4(t,in2,in3,in4,in5,in6)
		%CONSCOR_4
		%    OUT1 = CONSCOR_4(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:38
		ax__FRAME_5__dt_1_ = in6(7);
		ax__FRAME_20__dt_1_ = in6(20);
		ay__FRAME_5__dt_1_ = in6(8);
		ay__FRAME_20__dt_1_ = in6(21);
		az__FRAME_5__dt_1_ = in6(9);
		az__FRAME_20__dt_1_ = in6(22);
		rqw__FRAME_5__dt_0_ = in6(10);
		rqw__FRAME_20__dt_0_ = in6(23);
		rqx__FRAME_5__dt_0_ = in6(11);
		rqx__FRAME_20__dt_0_ = in6(24);
		rqy__FRAME_5__dt_0_ = in6(12);
		rqy__FRAME_20__dt_0_ = in6(25);
		rqz__FRAME_5__dt_0_ = in6(13);
		rqz__FRAME_20__dt_0_ = in6(26);
		t2 = (ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t3 = (ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t4 = (az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t5 = (ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t6 = (ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t7 = (ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t8 = (az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t9 = (ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t10 = (az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t11 = (ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t12 = (az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t13 = (ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t14 = (az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t15 = (ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t16 = (ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t17 = (ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t18 = (ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t19 = (az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t20 = t2+t3+t4;
		t41 = (az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t21 = t5+t6-t41;
		t22 = (rqx__FRAME_20__dt_0_.*t21)./2.0;
		t36 = (ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t23 = t7+t8-t36;
		t24 = (rqy__FRAME_20__dt_0_.*t23)./2.0;
		t34 = (ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t25 = t9+t10-t34;
		t26 = (rqz__FRAME_20__dt_0_.*t25)./2.0;
		t27 = t17+t18+t19;
		t47 = (az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t28 = t15+t16-t47;
		t29 = (rqx__FRAME_5__dt_0_.*t28)./2.0;
		t46 = (ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t30 = t13+t14-t46;
		t31 = (rqy__FRAME_5__dt_0_.*t30)./2.0;
		t45 = (ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t32 = t11+t12-t45;
		t33 = (rqz__FRAME_5__dt_0_.*t32)./2.0;
		t35 = (az__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t37 = (ax__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t38 = (ay__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t39 = (az__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t40 = (ay__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t42 = (ay__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t43 = (az__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t63 = (ax__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t44 = t2+t3+t4+t42+t43-t63;
		t48 = (az__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t49 = (ax__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t50 = (ay__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t51 = (az__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t52 = (ay__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t53 = (ay__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t54 = (az__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t66 = (ax__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t55 = t17+t18+t19+t53+t54-t66;
		t56 = (rqw__FRAME_20__dt_0_.*t25)./2.0;
		t57 = (rqx__FRAME_20__dt_0_.*t23)./2.0;
		t58 = (rqz__FRAME_20__dt_0_.*t20)./2.0;
		t59 = (rqy__FRAME_5__dt_0_.*t28)./2.0;
		t60 = t9+t10-t34+t35-(ax__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0-(ay__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t61 = t7+t8-t36+t37+t38+t39;
		t62 = t5+t6+t40-t41-(ax__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t64 = t11+t12-t45+t48-(ax__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0-(ay__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t65 = t13+t14-t46+t49+t50+t51;
		t67 = (rqw__FRAME_20__dt_0_.*t23)./2.0;
		t68 = (rqx__FRAME_5__dt_0_.*t32)./2.0;
		t69 = (rqy__FRAME_20__dt_0_.*t20)./2.0;
		t70 = (rqz__FRAME_20__dt_0_.*t21)./2.0;
		out1 = -ay__FRAME_5__dt_1_.*(-t56-t57-t58-t59+(rqw__FRAME_5__dt_0_.*t32)./2.0+(rqw__FRAME_5__dt_0_.*t64)./2.0+(rqx__FRAME_5__dt_0_.*t30)./2.0+(rqx__FRAME_5__dt_0_.*t65)./2.0+(rqy__FRAME_20__dt_0_.*t21)./2.0+(rqz__FRAME_5__dt_0_.*t27)./2.0+(rqz__FRAME_5__dt_0_.*t55)./2.0-(rqy__FRAME_5__dt_0_.*(t15+t16-t47+t52-(ax__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0))./2.0)+az__FRAME_5__dt_1_.*(-t67-t68-t69-t70+(rqw__FRAME_5__dt_0_.*t30)./2.0+(rqw__FRAME_5__dt_0_.*t65)./2.0-(rqx__FRAME_5__dt_0_.*t64)./2.0+(rqx__FRAME_20__dt_0_.*t25)./2.0+(rqy__FRAME_5__dt_0_.*t27)./2.0+(rqy__FRAME_5__dt_0_.*t55)./2.0+(rqz__FRAME_5__dt_0_.*t28)./2.0+(rqz__FRAME_5__dt_0_.*(t15+t16-t47+t52-(ax__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0))./2.0)+ay__FRAME_20__dt_1_.*(t56+t57+t58+t59-(rqw__FRAME_5__dt_0_.*t32)./2.0+(rqw__FRAME_20__dt_0_.*t60)./2.0-(rqx__FRAME_5__dt_0_.*t30)./2.0+(rqx__FRAME_20__dt_0_.*t61)./2.0-(rqy__FRAME_20__dt_0_.*t21)./2.0-(rqy__FRAME_20__dt_0_.*t62)./2.0-(rqz__FRAME_5__dt_0_.*t27)./2.0+(rqz__FRAME_20__dt_0_.*t44)./2.0)-az__FRAME_20__dt_1_.*(t67+t68+t69+t70-(rqw__FRAME_5__dt_0_.*t30)./2.0+(rqw__FRAME_20__dt_0_.*t61)./2.0-(rqx__FRAME_20__dt_0_.*t25)./2.0-(rqx__FRAME_20__dt_0_.*t60)./2.0-(rqy__FRAME_5__dt_0_.*t27)./2.0+(rqy__FRAME_20__dt_0_.*t44)./2.0-(rqz__FRAME_5__dt_0_.*t28)./2.0+(rqz__FRAME_20__dt_0_.*t62)./2.0)-ax__FRAME_20__dt_1_.*(t22+t24+t26+t29+t31+t33+(rqz__FRAME_20__dt_0_.*(t9+t10+t35-(ax__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0-(ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0-(ay__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0))./2.0+(rqx__FRAME_20__dt_0_.*(t5+t6+t40-(ax__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0))./2.0-(rqw__FRAME_5__dt_0_.*t27)./2.0-(rqw__FRAME_20__dt_0_.*t20)./2.0-(rqw__FRAME_20__dt_0_.*t44)./2.0+(rqy__FRAME_20__dt_0_.*(t7+t8+t37+t38+t39-(ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0))./2.0)+ax__FRAME_5__dt_1_.*(t22+t24+t26+t29+t31+t33+(rqz__FRAME_5__dt_0_.*(t11+t12+t48-(ax__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0-(ay__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0-(ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0))./2.0+(rqx__FRAME_5__dt_0_.*(t15+t16+t52-(ax__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0))./2.0-(rqw__FRAME_5__dt_0_.*t27)./2.0-(rqw__FRAME_5__dt_0_.*t55)./2.0-(rqw__FRAME_20__dt_0_.*t20)./2.0+(rqy__FRAME_5__dt_0_.*(t13+t14+t49+t50+t51-(ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0))./2.0);
		end
		function out1 = ConsCor_5(t,in2,in3,in4,in5,in6)
		%CONSCOR_5
		%    OUT1 = CONSCOR_5(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:41
		ax__FRAME_5__dt_1_ = in6(7);
		ax__FRAME_20__dt_1_ = in6(20);
		ay__FRAME_5__dt_1_ = in6(8);
		ay__FRAME_20__dt_1_ = in6(21);
		az__FRAME_5__dt_1_ = in6(9);
		az__FRAME_20__dt_1_ = in6(22);
		rqw__FRAME_5__dt_0_ = in6(10);
		rqw__FRAME_20__dt_0_ = in6(23);
		rqx__FRAME_5__dt_0_ = in6(11);
		rqx__FRAME_20__dt_0_ = in6(24);
		rqy__FRAME_5__dt_0_ = in6(12);
		rqy__FRAME_20__dt_0_ = in6(25);
		rqz__FRAME_5__dt_0_ = in6(13);
		rqz__FRAME_20__dt_0_ = in6(26);
		t2 = (ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t3 = (az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t4 = (ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t5 = (ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t6 = (ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t7 = (ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t8 = (ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t9 = (az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t10 = (ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t11 = (az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t12 = (ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t13 = (az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t14 = (ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t15 = (az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t16 = (ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t17 = (ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t18 = (ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t19 = (ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t20 = (az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t34 = (ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t21 = t2+t3-t34;
		t22 = (rqw__FRAME_20__dt_0_.*t21)./2.0;
		t23 = -t4+t10+t11;
		t24 = (rqx__FRAME_20__dt_0_.*t23)./2.0;
		t47 = (ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t25 = t12+t13-t47;
		t41 = (az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t26 = t5+t6-t41;
		t46 = (ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t27 = t14+t15-t46;
		t28 = t7+t8+t9;
		t29 = (rqz__FRAME_20__dt_0_.*t28)./2.0;
		t45 = (az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t30 = t16+t17-t45;
		t31 = (rqy__FRAME_5__dt_0_.*t30)./2.0;
		t32 = t18+t19+t20;
		t33 = (ax__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t35 = (ay__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t36 = (az__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t37 = (ax__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t38 = (ay__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t64 = (az__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t39 = t4-t10-t11+t37+t38-t64;
		t40 = (ax__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t42 = (ax__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t43 = (az__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t65 = (ay__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t44 = t7+t8+t9+t42+t43-t65;
		t48 = (ax__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t49 = (ay__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t50 = (az__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t51 = (az__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t52 = (ax__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t53 = (ax__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t54 = (az__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t69 = (ay__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t55 = t18+t19+t20+t53+t54-t69;
		t56 = (rqx__FRAME_20__dt_0_.*t26)./2.0;
		t57 = (rqy__FRAME_20__dt_0_.*t23)./2.0;
		t58 = (rqz__FRAME_20__dt_0_.*t21)./2.0;
		t59 = (rqx__FRAME_5__dt_0_.*t30)./2.0;
		t60 = (rqy__FRAME_5__dt_0_.*t27)./2.0;
		t61 = (rqz__FRAME_5__dt_0_.*t25)./2.0;
		t62 = t2+t3+t33-t34+t35+t36;
		t63 = t5+t6+t40-t41-(ay__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t66 = t12+t13-t47+t48+t49+t50;
		t67 = t16+t17-t45+t52-(ay__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t68 = t14+t15-t46+t51-(ax__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0-(ay__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t70 = (rqw__FRAME_20__dt_0_.*t26)./2.0;
		t71 = (rqy__FRAME_20__dt_0_.*t21)./2.0;
		t72 = (rqx__FRAME_20__dt_0_.*t28)./2.0;
		t73 = (rqz__FRAME_5__dt_0_.*t27)./2.0;
		out1 = az__FRAME_20__dt_1_.*(t70+t71+t72+t73-(rqw__FRAME_5__dt_0_.*t30)./2.0+(rqw__FRAME_20__dt_0_.*t63)./2.0-(rqx__FRAME_5__dt_0_.*t32)./2.0+(rqx__FRAME_20__dt_0_.*t44)./2.0-(rqy__FRAME_5__dt_0_.*t25)./2.0+(rqy__FRAME_20__dt_0_.*t62)./2.0-(rqz__FRAME_20__dt_0_.*t23)./2.0+(rqz__FRAME_20__dt_0_.*t39)./2.0)-ax__FRAME_20__dt_1_.*(t22+t24+t29+t31-(rqy__FRAME_20__dt_0_.*(t5+t6+t40-(ay__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0))./2.0-(rqw__FRAME_5__dt_0_.*t25)./2.0-(rqx__FRAME_5__dt_0_.*t27)./2.0-(rqx__FRAME_20__dt_0_.*t39)./2.0-(rqy__FRAME_20__dt_0_.*t26)./2.0-(rqz__FRAME_5__dt_0_.*t32)./2.0+(rqz__FRAME_20__dt_0_.*t44)./2.0+(rqw__FRAME_20__dt_0_.*(t2+t3+t33+t35+t36-(ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0))./2.0)+ay__FRAME_5__dt_1_.*(t56+t57+t58+t59+t60+t61-(rqw__FRAME_5__dt_0_.*t32)./2.0-(rqw__FRAME_5__dt_0_.*t55)./2.0-(rqw__FRAME_20__dt_0_.*t28)./2.0+(rqx__FRAME_5__dt_0_.*t67)./2.0+(rqy__FRAME_5__dt_0_.*t68)./2.0+(rqz__FRAME_5__dt_0_.*t66)./2.0)-ay__FRAME_20__dt_1_.*(t56+t57+t58+t59+t60+t61-(rqw__FRAME_5__dt_0_.*t32)./2.0-(rqw__FRAME_20__dt_0_.*t28)./2.0-(rqw__FRAME_20__dt_0_.*t44)./2.0+(rqx__FRAME_20__dt_0_.*t63)./2.0-(rqy__FRAME_20__dt_0_.*t39)./2.0+(rqz__FRAME_20__dt_0_.*t62)./2.0)-az__FRAME_5__dt_1_.*(-t70-t71-t72-t73+(rqw__FRAME_5__dt_0_.*t30)./2.0+(rqw__FRAME_5__dt_0_.*t67)./2.0+(rqx__FRAME_5__dt_0_.*t32)./2.0+(rqx__FRAME_5__dt_0_.*t55)./2.0+(rqy__FRAME_5__dt_0_.*t25)./2.0+(rqy__FRAME_5__dt_0_.*t66)./2.0-(rqz__FRAME_5__dt_0_.*t68)./2.0+(rqz__FRAME_20__dt_0_.*t23)./2.0)+ax__FRAME_5__dt_1_.*(-t22-t24-t29-t31+(rqx__FRAME_5__dt_0_.*(t14+t15+t51-(ax__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0-(ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0-(ay__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0))./2.0-(rqy__FRAME_5__dt_0_.*(t16+t17+t52-(ay__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0))./2.0+(rqw__FRAME_5__dt_0_.*t25)./2.0+(rqx__FRAME_5__dt_0_.*t27)./2.0+(rqy__FRAME_20__dt_0_.*t26)./2.0+(rqz__FRAME_5__dt_0_.*t32)./2.0+(rqz__FRAME_5__dt_0_.*t55)./2.0+(rqw__FRAME_5__dt_0_.*(t12+t13+t48+t49+t50-(ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0))./2.0);
		end
		function out1 = ConsCor_6(t,in2,in3,in4,in5,in6)
		%CONSCOR_6
		%    OUT1 = CONSCOR_6(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:44
		ax__FRAME_5__dt_1_ = in6(7);
		ax__FRAME_20__dt_1_ = in6(20);
		ay__FRAME_5__dt_1_ = in6(8);
		ay__FRAME_20__dt_1_ = in6(21);
		az__FRAME_5__dt_1_ = in6(9);
		az__FRAME_20__dt_1_ = in6(22);
		rqw__FRAME_5__dt_0_ = in6(10);
		rqw__FRAME_20__dt_0_ = in6(23);
		rqx__FRAME_5__dt_0_ = in6(11);
		rqx__FRAME_20__dt_0_ = in6(24);
		rqy__FRAME_5__dt_0_ = in6(12);
		rqy__FRAME_20__dt_0_ = in6(25);
		rqz__FRAME_5__dt_0_ = in6(13);
		rqz__FRAME_20__dt_0_ = in6(26);
		t2 = (ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t3 = (ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t4 = (az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t5 = (ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t6 = (ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t7 = (az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t8 = (ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t9 = (ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t10 = (ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t11 = (az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t12 = (ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t13 = (az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t14 = (ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t15 = (ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t16 = (ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t17 = (ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t18 = (az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t19 = (ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t20 = (az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t21 = -t2+t19+t20;
		t22 = (rqw__FRAME_20__dt_0_.*t21)./2.0;
		t37 = (ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t23 = t3+t4-t37;
		t33 = (ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t24 = t12+t13-t33;
		t34 = (ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t25 = t10+t11-t34;
		t26 = (rqx__FRAME_5__dt_0_.*t25)./2.0;
		t27 = t5+t6+t7;
		t28 = (rqy__FRAME_20__dt_0_.*t27)./2.0;
		t40 = (az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t29 = t8+t9-t40;
		t30 = (rqz__FRAME_20__dt_0_.*t29)./2.0;
		t31 = t16+t17+t18;
		t35 = (az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t32 = t14+t15-t35;
		t36 = (ax__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t38 = (ax__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t39 = (ay__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t41 = (az__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t42 = (ax__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t43 = (az__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t44 = (ax__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
		t45 = (ay__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
		t65 = (az__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
		t46 = t5+t6+t7+t44+t45-t65;
		t47 = (ax__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t66 = (ay__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t67 = (az__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t48 = t10+t11-t34+t47-t66-t67;
		t49 = (ax__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t50 = (ay__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t51 = (az__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t52 = t14+t15-t35+t49+t50+t51;
		t53 = (ay__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
		t68 = (ax__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t69 = (az__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t54 = t12+t13-t33+t53-t68-t69;
		t55 = (ax__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
		t56 = (ay__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
		t70 = (az__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
		t57 = t16+t17+t18+t55+t56-t70;
		t58 = (rqw__FRAME_20__dt_0_.*t29)./2.0;
		t59 = (rqy__FRAME_20__dt_0_.*t23)./2.0;
		t60 = (rqx__FRAME_20__dt_0_.*t27)./2.0;
		t61 = (rqz__FRAME_5__dt_0_.*t24)./2.0;
		t62 = t3+t4+t36-t37-(ay__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t63 = t2-t19-t20+t42+t43-(ay__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
		t64 = t8+t9+t38+t39-t40+t41;
		t71 = (rqx__FRAME_20__dt_0_.*t29)./2.0;
		t72 = (rqy__FRAME_20__dt_0_.*t21)./2.0;
		t73 = (rqz__FRAME_20__dt_0_.*t23)./2.0;
		t74 = (rqx__FRAME_5__dt_0_.*t32)./2.0;
		t75 = (rqy__FRAME_5__dt_0_.*t24)./2.0;
		t76 = (rqz__FRAME_5__dt_0_.*t25)./2.0;
		out1 = -ay__FRAME_20__dt_1_.*(t58+t59+t60+t61-(rqw__FRAME_5__dt_0_.*t32)./2.0+(rqw__FRAME_20__dt_0_.*t64)./2.0-(rqx__FRAME_5__dt_0_.*t31)./2.0+(rqx__FRAME_20__dt_0_.*t46)./2.0-(rqy__FRAME_5__dt_0_.*t25)./2.0+(rqy__FRAME_20__dt_0_.*t62)./2.0-(rqz__FRAME_20__dt_0_.*t21)./2.0+(rqz__FRAME_20__dt_0_.*t63)./2.0)+az__FRAME_5__dt_1_.*(t71+t72+t73+t74+t75+t76+(rqy__FRAME_5__dt_0_.*(t12+t13-t33+t53-t68-t69))./2.0-(rqw__FRAME_5__dt_0_.*t31)./2.0-(rqw__FRAME_5__dt_0_.*t57)./2.0-(rqw__FRAME_20__dt_0_.*t27)./2.0+(rqx__FRAME_5__dt_0_.*t52)./2.0+(rqz__FRAME_5__dt_0_.*t48)./2.0)+ax__FRAME_20__dt_1_.*(t22+t26+t28+t30-(rqw__FRAME_20__dt_0_.*(t2+t42+t43-(ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0-(ay__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0))./2.0-(rqx__FRAME_20__dt_0_.*(t3+t4+t36-(ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0-(ay__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0))./2.0-(rqw__FRAME_5__dt_0_.*t24)./2.0-(rqx__FRAME_20__dt_0_.*t23)./2.0-(rqy__FRAME_5__dt_0_.*t31)./2.0+(rqy__FRAME_20__dt_0_.*t46)./2.0-(rqz__FRAME_5__dt_0_.*t32)./2.0+(rqz__FRAME_20__dt_0_.*(t8+t9+t38+t39+t41-(az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0))./2.0)-az__FRAME_20__dt_1_.*(t71+t72+t73+t74+t75+t76-(rqw__FRAME_5__dt_0_.*t31)./2.0-(rqw__FRAME_20__dt_0_.*t27)./2.0-(rqw__FRAME_20__dt_0_.*t46)./2.0+(rqx__FRAME_20__dt_0_.*t64)./2.0-(rqy__FRAME_20__dt_0_.*t63)./2.0+(rqz__FRAME_20__dt_0_.*t62)./2.0)-ax__FRAME_5__dt_1_.*(-t22-t26-t28-t30+(rqw__FRAME_5__dt_0_.*t24)./2.0+(rqw__FRAME_5__dt_0_.*t54)./2.0-(rqx__FRAME_5__dt_0_.*t48)./2.0+(rqx__FRAME_20__dt_0_.*t23)./2.0+(rqy__FRAME_5__dt_0_.*t31)./2.0+(rqy__FRAME_5__dt_0_.*t57)./2.0+(rqz__FRAME_5__dt_0_.*t32)./2.0+(rqz__FRAME_5__dt_0_.*t52)./2.0)+ay__FRAME_5__dt_1_.*(-t58-t59-t60-t61+(rqw__FRAME_5__dt_0_.*t32)./2.0+(rqw__FRAME_5__dt_0_.*t52)./2.0+(rqx__FRAME_5__dt_0_.*t31)./2.0+(rqx__FRAME_5__dt_0_.*t57)./2.0+(rqy__FRAME_5__dt_0_.*t25)./2.0+(rqy__FRAME_5__dt_0_.*t48)./2.0-(rqz__FRAME_5__dt_0_.*t54)./2.0+(rqz__FRAME_20__dt_0_.*t21)./2.0);
		end
		function out1 = ConsCor_7(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSCOR_7
		%    OUT1 = CONSCOR_7(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:20:19
		q1dev__dt_0_ = in2(1,:);
		q1dev__dt_1_ = in2(10,:);
		q1flex__dt_0_ = in2(2,:);
		q1flex__dt_1_ = in2(11,:);
		q1sup__dt_0_ = in2(3,:);
		q1sup__dt_1_ = in2(12,:);
		xx3__dt_0_ = in5(1,:);
		zz3__dt_0_ = in5(3,:);
		t2 = q1dev__dt_0_./2.0;
		t3 = q1flex__dt_0_./2.0;
		t4 = q1sup__dt_0_./2.0;
		t5 = cos(t2);
		t6 = cos(t3);
		t7 = sin(t4);
		t8 = cos(t4);
		t9 = sin(t2);
		t10 = sin(t3);
		t11 = t8.^2;
		t12 = t5.^2;
		t13 = t6.^2;
		t14 = xx3__dt_0_./2.0;
		t15 = zz3__dt_0_./2.0;
		t16 = cos(q1dev__dt_0_);
		t17 = cos(q1sup__dt_0_);
		t18 = t11.*t12.*2.0;
		t19 = t11.*t13.*2.0;
		t20 = t11.^2;
		t21 = t7.*t9.*t10.*2.0;
		t22 = t5.*t6.*t8.*t11.*2.0;
		t29 = t11.*3.0;
		t30 = t7.*t9.*t10.*t11.*2.0;
		t23 = -t12-t13+t18+t19+t20+t21+t22-t29-t30+2.0;
		t24 = 1.0./t23;
		t25 = cos(t14);
		t26 = sin(t15);
		t27 = sin(q1sup__dt_0_);
		t28 = t6.*t8.*t9.*t11;
		t31 = cos(t15);
		t32 = sin(t14);
		t33 = t5.*t9.*t11.*2.0;
		t34 = t5.*t7.*t10.*t11;
		t61 = t5.*t7.*t10;
		t126 = t5.*t9;
		t35 = t28+t33+t34-t61-t126;
		t36 = sin(q1flex__dt_0_);
		t37 = sin(q1dev__dt_0_);
		t38 = 1.0./t23.^2;
		t39 = t5.*t6.*t7.*(3.0./4.0);
		t40 = (t8.*t9.*t10.*t11)./2.0;
		t41 = (t5.*t6.*t7.*t17)./4.0;
		t42 = cos(q1flex__dt_0_);
		t43 = t7.*t9.*t10.*(3.0./4.0);
		t44 = (t5.*t6.*t8.*t11)./2.0;
		t45 = (t7.*t9.*t10.*t17)./4.0;
		t46 = t37./2.0;
		t47 = (t27.*t36)./2.0;
		t48 = t6.*t8.*t9.*2.0;
		t49 = t5.*t7.*t10.*(3.0./2.0);
		t50 = (t5.*t7.*t10.*t17)./2.0;
		t51 = -t28+t46+t47+t48+t49+t50;
		t52 = t5.*t8.*t10.*t11;
		t53 = t6.*t10.*t11.*2.0;
		t54 = t6.*t7.*t9.*t11;
		t62 = t6.*t7.*t9;
		t139 = t6.*t10;
		t55 = t52+t53+t54-t62-t139;
		t56 = t36./2.0;
		t57 = (t27.*t37)./2.0;
		t58 = t5.*t8.*t10.*2.0;
		t59 = t6.*t7.*t9.*(3.0./2.0);
		t60 = (t6.*t7.*t9.*t17)./2.0;
		t63 = -t52+t56+t57+t58+t59+t60;
		t64 = t7.*t8.*t11.*2.0;
		t65 = t7.*t8.*t12.*2.0;
		t66 = t7.*t8.*t13.*2.0;
		t67 = t8.*t9.*t10.*t11;
		t68 = t5.*t6.*t7.*t11.*3.0;
		t69 = t7.^2;
		t115 = t7.*t8.*3.0;
		t116 = t8.*t9.*t10;
		t117 = t8.*t9.*t10.*t69.*2.0;
		t70 = t64+t65+t66+t67+t68-t115-t116-t117;
		t71 = (t6.*t8.*t9)./4.0;
		t72 = (t5.*t7.*t10)./4.0;
		t73 = (t5.*t8.*t10)./4.0;
		t74 = (t6.*t7.*t9)./4.0;
		t75 = q1flex__dt_1_.*t17;
		t76 = (q1dev__dt_1_.*t27)./2.0;
		t77 = (q1sup__dt_1_.*t37)./2.0;
		t78 = (q1sup__dt_1_.*t27.*t36)./2.0;
		t79 = q1flex__dt_1_.*t5.*t6.*t8;
		t80 = (q1dev__dt_1_.*t5.*t6.*t7)./2.0;
		t81 = q1sup__dt_1_.*t6.*t8.*t9.*2.0;
		t82 = q1dev__dt_1_.*t8.*t9.*t10;
		t83 = q1sup__dt_1_.*t5.*t7.*t10.*(3.0./2.0);
		t84 = (q1dev__dt_1_.*t5.*t6.*t7.*t17)./2.0;
		t85 = (q1sup__dt_1_.*t5.*t7.*t10.*t17)./2.0;
		t86 = q1flex__dt_1_.*t7.*t9.*t10.*t17.*(3.0./2.0);
		t108 = q1flex__dt_1_.*t20.*2.0;
		t109 = (q1flex__dt_1_.*t16.*t17)./2.0;
		t110 = (q1flex__dt_1_.*t17.*t42)./2.0;
		t111 = (q1flex__dt_1_.*t7.*t9.*t10)./2.0;
		t112 = q1flex__dt_1_.*t5.*t6.*t8.*t11.*3.0;
		t113 = q1sup__dt_1_.*t6.*t8.*t9.*t11;
		t114 = q1dev__dt_1_.*t8.*t9.*t10.*t11;
		t87 = t75+t76+t77+t78+t79+t80+t81+t82+t83+t84+t85+t86-t108-t109-t110-t111-t112-t113-t114;
		t88 = q1dev__dt_1_.*t17;
		t89 = (q1flex__dt_1_.*t27)./2.0;
		t90 = (q1sup__dt_1_.*t36)./2.0;
		t91 = (q1sup__dt_1_.*t27.*t37)./2.0;
		t92 = q1dev__dt_1_.*t5.*t6.*t8;
		t93 = (q1flex__dt_1_.*t5.*t6.*t7)./2.0;
		t94 = q1sup__dt_1_.*t5.*t8.*t10.*2.0;
		t95 = q1flex__dt_1_.*t8.*t9.*t10;
		t96 = q1sup__dt_1_.*t6.*t7.*t9.*(3.0./2.0);
		t97 = (q1flex__dt_1_.*t5.*t6.*t7.*t17)./2.0;
		t98 = (q1sup__dt_1_.*t6.*t7.*t9.*t17)./2.0;
		t99 = q1dev__dt_1_.*t7.*t9.*t10.*t17.*(3.0./2.0);
		t101 = q1dev__dt_1_.*t20.*2.0;
		t102 = (q1dev__dt_1_.*t16.*t17)./2.0;
		t103 = (q1dev__dt_1_.*t17.*t42)./2.0;
		t104 = (q1dev__dt_1_.*t7.*t9.*t10)./2.0;
		t105 = q1dev__dt_1_.*t5.*t6.*t8.*t11.*3.0;
		t106 = q1sup__dt_1_.*t5.*t8.*t10.*t11;
		t107 = q1flex__dt_1_.*t8.*t9.*t10.*t11;
		t100 = t88+t89+t90+t91+t92+t93+t94+t95+t96+t97+t98+t99-t101-t102-t103-t104-t105-t106-t107;
		t118 = (t5.*t6.*t7)./4.0;
		t119 = (t8.*t9.*t10)./4.0;
		t120 = t118+t119;
		t121 = t71+t72;
		t122 = (t16.*t27)./2.0;
		t123 = (t27.*t42)./2.0;
		t124 = (t5.*t6.*t8)./4.0;
		t125 = (t17.*t37)./2.0;
		t127 = (t16.*t17)./2.0;
		t128 = (t17.*t42)./2.0;
		t129 = t20.*2.0;
		t130 = (t7.*t9.*t10)./2.0;
		t131 = t5.*t6.*t8.*t11.*3.0;
		t140 = t5.*t6.*t8;
		t141 = t7.*t9.*t10.*t17.*(3.0./2.0);
		t132 = -t17+t127+t128+t129+t130+t131-t140-t141;
		t133 = t27./2.0;
		t134 = (t5.*t6.*t7)./2.0;
		t135 = (t5.*t6.*t7.*t17)./2.0;
		t136 = -t67+t116+t133+t134+t135;
		t137 = (t7.*t9.*t10)./4.0;
		t138 = (t17.*t36)./2.0;
		t142 = (t24.*t25.*t31.*t87)./8.0;
		t175 = (t24.*t26.*t32.*t100)./8.0;
		t143 = t142-t175;
		t144 = (t24.*t26.*t32.*t87)./8.0;
		t174 = (t24.*t25.*t31.*t100)./8.0;
		t145 = t144-t174;
		t146 = t124+t137;
		t147 = t73+t74;
		t148 = t7.*t8.*t11.*4.0;
		t149 = t5.*t6.*t7.*t11.*(9.0./2.0);
		t150 = t8.*t9.*t10.*t17.*(3.0./4.0);
		t151 = -t27-t119+t122+t123-t134+t148+t149+t150-t7.*t9.*t10.*t27.*(3.0./2.0);
		t152 = t17./2.0;
		t153 = t7.*t9.*t10.*t11.*(3.0./2.0);
		t154 = (t5.*t6.*t8.*t17)./4.0;
		t155 = t124-t130+t152+t153+t154-(t5.*t6.*t7.*t27)./2.0;
		t156 = (t5.*t8.*t10.*t11)./2.0;
		t157 = (t6.*t7.*t9.*t17)./4.0;
		t165 = (t5.*t8.*t10)./2.0;
		t158 = t74+t156+t157-t165;
		t159 = t6.*t8.*t9.*t11.*(3.0./2.0);
		t160 = t5.*t7.*t10.*t17.*(3.0./4.0);
		t162 = (t6.*t8.*t9)./2.0;
		t161 = -t72+t125+t159+t160-t162;
		t163 = (t6.*t8.*t9.*t11)./2.0;
		t164 = (t5.*t7.*t10.*t17)./4.0;
		t166 = t5.*t8.*t10.*t11.*(3.0./2.0);
		t167 = t6.*t7.*t9.*t17.*(3.0./4.0);
		t168 = (t24.*t25.*t31.*t132)./8.0;
		t169 = (t24.*t26.*t32.*t136)./8.0;
		t170 = t168+t169;
		t171 = (t24.*t25.*t31.*t136)./8.0;
		t172 = (t24.*t26.*t32.*t132)./8.0;
		t173 = t171+t172;
		t176 = (q1flex__dt_1_.*t5.*t8.*t10)./2.0;
		t177 = q1sup__dt_1_.*t5.*t6.*t7.*(3.0./4.0);
		t178 = (q1sup__dt_1_.*t8.*t9.*t10.*t11)./2.0;
		t179 = (q1sup__dt_1_.*t5.*t6.*t7.*t17)./4.0;
		t180 = (q1flex__dt_1_.*t6.*t8.*t9)./2.0;
		t181 = (q1dev__dt_1_.*t6.*t7.*t9)./4.0;
		t182 = (q1flex__dt_1_.*t5.*t7.*t10)./4.0;
		t183 = q1sup__dt_1_.*t7.*t9.*t10.*(3.0./4.0);
		t184 = (q1sup__dt_1_.*t5.*t6.*t8.*t11)./2.0;
		t185 = (q1sup__dt_1_.*t7.*t9.*t10.*t17)./4.0;
		out1 = q1flex__dt_1_.*(-q1dev__dt_1_.*t146+q1flex__dt_1_.*t120+q1sup__dt_1_.*t147-q1dev__dt_1_.*(t124+t137+(t24.*t25.*t26.*t161)./4.0-(t24.*t31.*t32.*t158)./4.0-(t25.*t26.*t35.*t38.*t132)./4.0+(t31.*t32.*t35.*t38.*t136)./4.0)+q1sup__dt_1_.*(t73+t74-(t24.*t25.*t26.*t151)./4.0-(t24.*t31.*t32.*t155)./4.0+(t25.*t26.*t38.*t70.*t132)./4.0-(t31.*t32.*t38.*t70.*t136)./4.0)+q1flex__dt_1_.*(t118+t119-(t24.*t25.*t26.*(-t74+t138-t165+t166+t167))./4.0+(t24.*t31.*t32.*(t72-t162+t163+t164))./4.0+(t25.*t26.*t38.*t55.*t132)./4.0-(t31.*t32.*t38.*t55.*t136)./4.0)+t24.*t87.*t173-t24.*t100.*t170+t24.*t132.*t145+t24.*t136.*t143+(t24.*t31.*t32.*(-t180+t181+t182+t183+t184+t185-(q1sup__dt_1_.*t42)./2.0-(q1dev__dt_1_.*t17.*t36)./2.0+(q1dev__dt_1_.*t5.*t8.*t10)./2.0-q1sup__dt_1_.*t5.*t6.*t8-q1dev__dt_1_.*t5.*t8.*t10.*t11.*(3.0./2.0)-q1dev__dt_1_.*t6.*t7.*t9.*t17.*(3.0./4.0)+(q1flex__dt_1_.*t6.*t8.*t9.*t11)./2.0+(q1flex__dt_1_.*t5.*t7.*t10.*t17)./4.0))./4.0-(t24.*t25.*t26.*(-t176+t177+t178+t179+(q1flex__dt_1_.*t17.*t36)./2.0+(q1sup__dt_1_.*t27.*t42)./2.0-(q1dev__dt_1_.*t5.*t7.*t10)./4.0+(q1dev__dt_1_.*t6.*t8.*t9)./2.0-(q1flex__dt_1_.*t6.*t7.*t9)./4.0-q1sup__dt_1_.*t8.*t9.*t10-(q1dev__dt_1_.*t6.*t8.*t9.*t11)./2.0-(q1dev__dt_1_.*t5.*t7.*t10.*t17)./4.0+q1flex__dt_1_.*t5.*t8.*t10.*t11.*(3.0./2.0)+q1flex__dt_1_.*t6.*t7.*t9.*t17.*(3.0./4.0)))./4.0-(t25.*t26.*t38.*t55.*t87)./4.0-(t31.*t32.*t38.*t55.*t100)./4.0)+q1sup__dt_1_.*(q1flex__dt_1_.*(t73+t74+(t24.*t31.*t32.*(t42.*(-1.0./2.0)+t43+t44+t45-t5.*t6.*t8))./4.0-(t24.*t25.*t26.*(t39+t40+t41+t123-t8.*t9.*t10))./4.0-(t25.*t26.*t38.*t51.*t55)./4.0-(t31.*t32.*t38.*t55.*t63)./4.0)+q1dev__dt_1_.*t121+q1flex__dt_1_.*t147+q1sup__dt_1_.*t120-q1sup__dt_1_.*(t5.*t6.*t7.*(-1.0./4.0)-(t8.*t9.*t10)./4.0+(t24.*t31.*t32.*(-t61+t125+t6.*t8.*t9.*(3.0./4.0)+t5.*t7.*t10.*t11.*(3.0./2.0)+(t6.*t8.*t9.*t17)./4.0-(t6.*t7.*t9.*t27)./2.0))./4.0+(t24.*t25.*t26.*(-t62+t138+t5.*t8.*t10.*(3.0./4.0)+t6.*t7.*t9.*t11.*(3.0./2.0)+(t5.*t8.*t10.*t17)./4.0-(t5.*t7.*t10.*t27)./2.0))./4.0+(t25.*t26.*t38.*t51.*t70)./4.0+(t31.*t32.*t38.*t63.*t70)./4.0)+q1dev__dt_1_.*(t71+t72+(t24.*t25.*t26.*(t16.*(-1.0./2.0)+t43+t44+t45-t5.*t6.*t8))./4.0-(t24.*t31.*t32.*(t39+t40+t41+t122-t8.*t9.*t10))./4.0-(t31.*t32.*t35.*t38.*(t56+t57+t58+t59+t60-t5.*t8.*t10.*t11))./4.0-(t25.*t26.*t35.*t38.*t51)./4.0)-t24.*t87.*((t24.*t26.*t32.*t51)./8.0-(t24.*t25.*t31.*t63)./8.0)+t24.*t100.*((t24.*t25.*t31.*t51)./8.0-(t24.*t26.*t32.*t63)./8.0)-t24.*t51.*t145+t24.*t63.*t143-(t24.*t31.*t32.*(-t80-t111-q1dev__dt_1_.*t27+(q1flex__dt_1_.*t17)./2.0+(q1dev__dt_1_.*t16.*t27)./2.0+(q1dev__dt_1_.*t27.*t42)./2.0+(q1sup__dt_1_.*t17.*t37)./2.0+q1dev__dt_1_.*t7.*t8.*t11.*4.0-(q1dev__dt_1_.*t8.*t9.*t10)./4.0+(q1flex__dt_1_.*t5.*t6.*t8)./4.0-q1sup__dt_1_.*t5.*t7.*t10+q1sup__dt_1_.*t6.*t8.*t9.*(3.0./4.0)+q1dev__dt_1_.*t5.*t6.*t7.*t11.*(9.0./2.0)+q1dev__dt_1_.*t8.*t9.*t10.*t17.*(3.0./4.0)-q1dev__dt_1_.*t7.*t9.*t10.*t27.*(3.0./2.0)+(q1flex__dt_1_.*t5.*t6.*t8.*t17)./4.0+q1flex__dt_1_.*t7.*t9.*t10.*t11.*(3.0./2.0)-(q1flex__dt_1_.*t5.*t6.*t7.*t27)./2.0+q1sup__dt_1_.*t5.*t7.*t10.*t11.*(3.0./2.0)+(q1sup__dt_1_.*t6.*t8.*t9.*t17)./4.0-(q1sup__dt_1_.*t6.*t7.*t9.*t27)./2.0))./4.0-(t24.*t25.*t26.*(-t93-t104+(q1dev__dt_1_.*t17)./2.0-q1flex__dt_1_.*t27+(q1flex__dt_1_.*t16.*t27)./2.0+(q1flex__dt_1_.*t27.*t42)./2.0+(q1sup__dt_1_.*t17.*t36)./2.0+(q1dev__dt_1_.*t5.*t6.*t8)./4.0+q1flex__dt_1_.*t7.*t8.*t11.*4.0-(q1flex__dt_1_.*t8.*t9.*t10)./4.0-q1sup__dt_1_.*t6.*t7.*t9+q1sup__dt_1_.*t5.*t8.*t10.*(3.0./4.0)+(q1dev__dt_1_.*t5.*t6.*t8.*t17)./4.0+q1dev__dt_1_.*t7.*t9.*t10.*t11.*(3.0./2.0)-(q1dev__dt_1_.*t5.*t6.*t7.*t27)./2.0+q1flex__dt_1_.*t5.*t6.*t7.*t11.*(9.0./2.0)+q1flex__dt_1_.*t8.*t9.*t10.*t17.*(3.0./4.0)-q1flex__dt_1_.*t7.*t9.*t10.*t27.*(3.0./2.0)+q1sup__dt_1_.*t6.*t7.*t9.*t11.*(3.0./2.0)+(q1sup__dt_1_.*t5.*t8.*t10.*t17)./4.0-(q1sup__dt_1_.*t5.*t7.*t10.*t27)./2.0))./4.0-(t25.*t26.*t38.*t70.*t87)./4.0-(t31.*t32.*t38.*t70.*t100)./4.0)-q1dev__dt_1_.*(-q1dev__dt_1_.*t120+q1flex__dt_1_.*t146-q1sup__dt_1_.*t121+q1flex__dt_1_.*(t124+t137+(t24.*t31.*t32.*(-t74+t138+t166+t167-(t5.*t8.*t10)./2.0))./4.0-(t24.*t25.*t26.*(t72+t163+t164-(t6.*t8.*t9)./2.0))./4.0+(t25.*t26.*t38.*t55.*t136)./4.0-(t31.*t32.*t38.*t55.*t132)./4.0)-q1dev__dt_1_.*(t118+t119+(t24.*t25.*t26.*t158)./4.0-(t24.*t31.*t32.*t161)./4.0-(t25.*t26.*t35.*t38.*t136)./4.0+(t31.*t32.*t35.*t38.*t132)./4.0)-q1sup__dt_1_.*(t71+t72-(t24.*t25.*t26.*t155)./4.0-(t24.*t31.*t32.*t151)./4.0-(t25.*t26.*t38.*t70.*t136)./4.0+(t31.*t32.*t38.*t70.*t132)./4.0)+t24.*t87.*t170-t24.*t100.*t173+t24.*t132.*t143+t24.*t136.*t145+(t24.*t31.*t32.*(t176+t177+t178+t179+(q1dev__dt_1_.*t17.*t37)./2.0+(q1sup__dt_1_.*t16.*t27)./2.0-(q1dev__dt_1_.*t5.*t7.*t10)./4.0-(q1dev__dt_1_.*t6.*t8.*t9)./2.0-(q1flex__dt_1_.*t6.*t7.*t9)./4.0-q1sup__dt_1_.*t8.*t9.*t10+q1dev__dt_1_.*t6.*t8.*t9.*t11.*(3.0./2.0)+q1dev__dt_1_.*t5.*t7.*t10.*t17.*(3.0./4.0)-(q1flex__dt_1_.*t5.*t8.*t10.*t11)./2.0-(q1flex__dt_1_.*t6.*t7.*t9.*t17)./4.0))./4.0-(t24.*t25.*t26.*(t180+t181+t182+t183+t184+t185-(q1sup__dt_1_.*t16)./2.0-(q1flex__dt_1_.*t17.*t37)./2.0-(q1dev__dt_1_.*t5.*t8.*t10)./2.0-q1sup__dt_1_.*t5.*t6.*t8+(q1dev__dt_1_.*t5.*t8.*t10.*t11)./2.0+(q1dev__dt_1_.*t6.*t7.*t9.*t17)./4.0-q1flex__dt_1_.*t6.*t8.*t9.*t11.*(3.0./2.0)-q1flex__dt_1_.*t5.*t7.*t10.*t17.*(3.0./4.0)))./4.0+(t25.*t26.*t35.*t38.*t87)./4.0+(t31.*t32.*t35.*t38.*t100)./4.0);
		end
		function out1 = ConsGF_1(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSGF_1
		%    OUT1 = CONSGF_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:57
		c = in3(3,:);
		l21 = in3(19,:);
		l22 = in3(20,:);
		l23 = in3(21,:);
		l24 = in3(22,:);
		l25 = in3(23,:);
		l11x = in3(7,:);
		l12x = in3(10,:);
		l12y = in3(11,:);
		l13x = in3(13,:);
		l12z = in3(12,:);
		l13y = in3(14,:);
		l13z = in3(15,:);
		l21x = in3(16,:);
		q1dev__dt_0_ = in2(1,:);
		q1dev__dt_1_ = in2(10,:);
		q1flex__dt_0_ = in2(2,:);
		q1flex__dt_1_ = in2(11,:);
		q1sup__dt_0_ = in2(3,:);
		q1sup__dt_1_ = in2(12,:);
		q2w1__dt_0_ = in2(4,:);
		q2w1__dt_1_ = in2(13,:);
		q2w2__dt_0_ = in2(5,:);
		q2w2__dt_1_ = in2(14,:);
		q2w3__dt_0_ = in2(6,:);
		q2w3__dt_1_ = in2(15,:);
		q2w4__dt_0_ = in2(7,:);
		q2w4__dt_1_ = in2(16,:);
		q2w5__dt_0_ = in2(8,:);
		q2w5__dt_1_ = in2(17,:);
		q2w6__dt_0_ = in2(9,:);
		q2w6__dt_1_ = in2(18,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf2Medz = in3(35,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		rf1Devz = in3(29,:);
		riw1 = in3(36,:);
		riw2 = in3(37,:);
		riw3 = in3(38,:);
		riw4 = in3(39,:);
		riw5 = in3(40,:);
		riw6 = in3(41,:);
		t2 = cos(q1sup__dt_0_);
		t3 = cos(rf1Devx);
		t4 = sin(rf1Devz);
		t5 = cos(rf1Devz);
		t6 = sin(rf1Devx);
		t7 = sin(rf1Devy);
		t8 = cos(q1flex__dt_0_);
		t9 = sin(q1flex__dt_0_);
		t10 = q1flex__dt_0_./2.0;
		t14 = q1dev__dt_0_./2.0;
		t11 = cos(t14);
		t12 = cos(t10);
		t13 = sin(t10);
		t15 = q1sup__dt_0_./2.0;
		t16 = cos(t15);
		t17 = sin(t14);
		t18 = sin(t15);
		t19 = sin(q1dev__dt_0_);
		t20 = cos(q1dev__dt_0_);
		t21 = sin(q1sup__dt_0_);
		t22 = q2w5__dt_0_+riw5;
		t23 = q2w6__dt_0_+riw6;
		t24 = q2w4__dt_0_+riw4;
		t25 = cos(t23);
		t26 = sin(t23);
		t27 = q2w3__dt_0_+riw3;
		t28 = cos(t24);
		t29 = cos(t22);
		t30 = l25+l13z;
		t31 = t29.*t30;
		t32 = sin(t22);
		t33 = l13x.*t25;
		t42 = l13y.*t26;
		t34 = l24+t33-t42;
		t43 = t32.*t34;
		t35 = t31-t43;
		t36 = sin(t24);
		t37 = l13y.*t25;
		t38 = l13x.*t26;
		t39 = l23+t37+t38;
		t40 = q2w2__dt_0_+riw2;
		t41 = sin(t27);
		t44 = t28.*t39;
		t56 = t35.*t36;
		t45 = l22+t44-t56;
		t46 = cos(t27);
		t47 = t28.*t35;
		t48 = t36.*t39;
		t49 = t47+t48;
		t50 = cos(rf2Medx);
		t51 = sin(rf2Medz);
		t52 = cos(rf2Medz);
		t53 = sin(rf2Medx);
		t54 = sin(rf2Medy);
		t55 = cos(t40);
		t57 = t41.*t45;
		t58 = t46.*t49;
		t59 = t57+t58;
		t60 = sin(t40);
		t61 = t41.*t49;
		t65 = t45.*t46;
		t62 = t61-t65;
		t63 = q2w1__dt_0_+riw1;
		t64 = t55.*t59;
		t66 = t64-t60.*t62;
		t67 = cos(t63);
		t68 = t33-t42;
		t69 = t37+t38;
		t70 = t36.*t68;
		t71 = t28.*t32.*t69;
		t72 = t70+t71;
		t73 = t28.*t68;
		t82 = t32.*t36.*t69;
		t74 = t73-t82;
		t75 = sin(t63);
		t76 = t50.*t51;
		t94 = t52.*t53.*t54;
		t77 = t76-t94;
		t78 = t51.*t53;
		t79 = t50.*t52.*t54;
		t80 = t78+t79;
		t81 = t46.*t72;
		t83 = t41.*t74;
		t84 = t81+t83;
		t85 = t41.*t72;
		t88 = t46.*t74;
		t86 = t85-t88;
		t87 = cos(rf2Medy);
		t89 = t55.*t86;
		t90 = t60.*t84;
		t91 = t89+t90;
		t92 = t59.*t60;
		t93 = l21-t61+t65;
		t95 = t60.*t93;
		t96 = t64+t95;
		t97 = t30.*t32;
		t98 = t29.*t34;
		t99 = t97+t98;
		t100 = t28.*t46.*t99;
		t106 = t36.*t41.*t99;
		t101 = t100-t106;
		t102 = t36.*t46.*t99;
		t103 = t28.*t41.*t99;
		t104 = t102+t103;
		t105 = t55.*t104;
		t107 = t60.*t101;
		t108 = t105+t107;
		t109 = t44-t56;
		t114 = t46.*t109;
		t110 = t61-t114;
		t111 = t41.*(t44-t56);
		t112 = t58+t111;
		t113 = t55.*t112;
		t115 = t113-t60.*t110;
		t117 = t55.*t93;
		t116 = t92-t117;
		t118 = t3.*t4;
		t129 = t5.*t6.*t7;
		t119 = t118-t129;
		t120 = t9.*t19;
		t121 = t8.*t20.*t21;
		t122 = t120+t121;
		t123 = t12.^2;
		t124 = t11.^2;
		t125 = cos(rf1Devy);
		t126 = t8.*t20;
		t127 = t9.*t19.*t21;
		t128 = t126+t127;
		t130 = t4.*t6;
		t131 = t3.*t5.*t7;
		t132 = t130+t131;
		t133 = t75.*t99;
		t134 = t133-t67.*t116;
		t135 = t9.*t20;
		t136 = t135-t8.*t19.*t21;
		t137 = l12x.*t2.*t19;
		t138 = t75.*t116;
		t139 = t67.*t99;
		t140 = t138+t139;
		t141 = t8.*t19;
		t142 = t141-t9.*t20.*t21;
		t143 = l12z.*t122;
		t144 = l12x.*t2.*t20;
		out1 = c.*(l11x-l21x-t119.*(t137+l12y.*(t123.*-2.0-t124.*2.0+t123.*t124.*4.0+t11.*t12.*t13.*t16.*t17.*t18.*8.0+1.0)-l12z.*t136)-t80.*t96+t77.*t134+t132.*(-l12x.*t21+l12y.*t2.*t9+l12z.*t2.*t8)-t52.*t87.*t140+t5.*t125.*(t143+t144-l12y.*t142))-c.*(q2w3__dt_1_.*(-t80.*(t92+t55.*t62)+t66.*t67.*t77+t52.*t66.*t75.*t87)-q2w1__dt_1_.*(t77.*t140+t52.*t87.*t134)-q1flex__dt_1_.*(t119.*(l12z.*t128-l12y.*(t12.*t13.*2.0-t12.*t13.*t124.*4.0-t11.*t13.^2.*t16.*t17.*t18.*4.0+t11.*t16.*t17.*t18.*t123.*4.0))+t132.*(l12y.*t2.*t8-l12z.*t2.*t9)+t5.*t125.*(l12y.*t122+l12z.*t142))+q2w2__dt_1_.*(-t80.*t116+t67.*t77.*t96+t52.*t75.*t87.*t96)+q2w4__dt_1_.*(-t80.*(t55.*t110+t60.*(t58+t41.*t109))+t67.*t77.*t115+t52.*t75.*t87.*t115)+q1dev__dt_1_.*(t119.*(t143+t144+l12y.*(t11.*t17.*2.0-t11.*t17.*t123.*4.0-t12.*t13.*t16.*t17.^2.*t18.*4.0+t12.*t13.*t16.*t18.*t124.*4.0))+t5.*t125.*(t137+l12y.*t128-l12z.*t136))+q2w6__dt_1_.*(t77.*(t67.*t91+t29.*t69.*t75)+t80.*(t55.*t84-t60.*t86)+t52.*t87.*(t75.*t91-t29.*t67.*t69))-q2w5__dt_1_.*(t77.*(t35.*t75+t67.*t108)+t80.*(t55.*t101-t60.*t104)-t52.*t87.*(t35.*t67-t75.*t108))+q1sup__dt_1_.*(t119.*(l12y.*(t11.*t12.*t13.*t16.^2.*t17.*4.0-t11.*t12.*t13.*t17.*t18.^2.*4.0)-l12x.*t19.*t21+l12z.*t2.*t8.*t19)+t132.*(l12x.*t2+l12y.*t9.*t21+l12z.*t8.*t21)-t5.*t125.*(-l12x.*t20.*t21+l12y.*t2.*t9.*t20+l12z.*t2.*t8.*t20)));
		end
		function out1 = ConsGF_2(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSGF_2
		%    OUT1 = CONSGF_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:18:54
		c = in3(3,:);
		l21 = in3(19,:);
		l22 = in3(20,:);
		l23 = in3(21,:);
		l24 = in3(22,:);
		l25 = in3(23,:);
		l11y = in3(8,:);
		l12x = in3(10,:);
		l12y = in3(11,:);
		l13x = in3(13,:);
		l12z = in3(12,:);
		l13y = in3(14,:);
		l13z = in3(15,:);
		l21y = in3(17,:);
		q1dev__dt_0_ = in2(1,:);
		q1dev__dt_1_ = in2(10,:);
		q1flex__dt_0_ = in2(2,:);
		q1flex__dt_1_ = in2(11,:);
		q1sup__dt_0_ = in2(3,:);
		q1sup__dt_1_ = in2(12,:);
		q2w1__dt_0_ = in2(4,:);
		q2w1__dt_1_ = in2(13,:);
		q2w2__dt_0_ = in2(5,:);
		q2w2__dt_1_ = in2(14,:);
		q2w3__dt_0_ = in2(6,:);
		q2w3__dt_1_ = in2(15,:);
		q2w4__dt_0_ = in2(7,:);
		q2w4__dt_1_ = in2(16,:);
		q2w5__dt_0_ = in2(8,:);
		q2w5__dt_1_ = in2(17,:);
		q2w6__dt_0_ = in2(9,:);
		q2w6__dt_1_ = in2(18,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf2Medz = in3(35,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		rf1Devz = in3(29,:);
		riw1 = in3(36,:);
		riw2 = in3(37,:);
		riw3 = in3(38,:);
		riw4 = in3(39,:);
		riw5 = in3(40,:);
		riw6 = in3(41,:);
		t2 = q2w5__dt_0_+riw5;
		t3 = q2w6__dt_0_+riw6;
		t4 = q2w4__dt_0_+riw4;
		t5 = cos(t3);
		t6 = sin(t3);
		t7 = q2w3__dt_0_+riw3;
		t8 = cos(t4);
		t9 = cos(t2);
		t10 = l25+l13z;
		t11 = t9.*t10;
		t12 = sin(t2);
		t13 = l13x.*t5;
		t22 = l13y.*t6;
		t14 = l24+t13-t22;
		t23 = t12.*t14;
		t15 = t11-t23;
		t16 = sin(t4);
		t17 = l13y.*t5;
		t18 = l13x.*t6;
		t19 = l23+t17+t18;
		t20 = q2w2__dt_0_+riw2;
		t21 = cos(t7);
		t24 = t8.*t19;
		t31 = t15.*t16;
		t25 = l22+t24-t31;
		t26 = sin(t7);
		t27 = t8.*t15;
		t28 = t16.*t19;
		t29 = t27+t28;
		t30 = sin(t20);
		t32 = t25.*t26;
		t33 = t21.*t29;
		t34 = t32+t33;
		t35 = cos(t20);
		t36 = t21.*t25;
		t68 = t26.*t29;
		t37 = l21+t36-t68;
		t38 = q2w1__dt_0_+riw1;
		t41 = rf2Medx./2.0;
		t39 = cos(t41);
		t43 = rf2Medz./2.0;
		t40 = cos(t43);
		t42 = t39.^2;
		t44 = t40.^2;
		t45 = rf2Medy./2.0;
		t48 = q1dev__dt_0_./2.0;
		t46 = cos(t48);
		t50 = q1flex__dt_0_./2.0;
		t47 = cos(t50);
		t49 = t46.^2;
		t51 = t47.^2;
		t52 = q1sup__dt_0_./2.0;
		t53 = sin(q1dev__dt_0_);
		t56 = rf1Devx./2.0;
		t54 = cos(t56);
		t58 = rf1Devz./2.0;
		t55 = cos(t58);
		t57 = t54.^2;
		t59 = t55.^2;
		t60 = rf1Devy./2.0;
		t61 = sin(q1sup__dt_0_);
		t62 = cos(q1flex__dt_0_);
		t63 = cos(q1sup__dt_0_);
		t64 = sin(q1flex__dt_0_);
		t65 = sin(rf2Medz);
		t66 = sin(t38);
		t67 = t30.*t34;
		t147 = t35.*t37;
		t69 = t67-t147;
		t70 = cos(t38);
		t71 = t10.*t12;
		t72 = t9.*t14;
		t73 = t71+t72;
		t74 = sin(rf1Devz);
		t75 = cos(q1dev__dt_0_);
		t76 = t53.*t64;
		t77 = t61.*t62.*t75;
		t78 = t76+t77;
		t79 = l12z.*t78;
		t80 = sin(t48);
		t81 = cos(t52);
		t82 = sin(t50);
		t83 = sin(t52);
		t84 = l12x.*t63.*t75;
		t85 = t57.*t59.*4.0;
		t86 = cos(t60);
		t87 = sin(t56);
		t88 = sin(t60);
		t89 = sin(t58);
		t90 = t54.*t55.*t86.*t87.*t88.*t89.*8.0;
		t96 = t57.*2.0;
		t97 = t59.*2.0;
		t91 = t85+t90-t96-t97+1.0;
		t92 = cos(rf1Devy);
		t93 = t64.*t75;
		t94 = t93-t53.*t61.*t62;
		t95 = l12x.*t53.*t63;
		t98 = cos(rf1Devz);
		t99 = sin(rf1Devx);
		t100 = t98.*t99;
		t101 = cos(rf1Devx);
		t102 = sin(rf1Devy);
		t104 = t74.*t101.*t102;
		t103 = t100-t104;
		t105 = t62.*t75;
		t106 = t53.*t61.*t64;
		t107 = t105+t106;
		t108 = t53.*t62;
		t109 = t108-t61.*t64.*t75;
		t110 = cos(rf2Medz);
		t111 = sin(rf2Medx);
		t112 = t110.*t111;
		t113 = cos(rf2Medx);
		t114 = sin(rf2Medy);
		t137 = t65.*t113.*t114;
		t115 = t112-t137;
		t116 = t34.*t35;
		t117 = t36-t68;
		t118 = t42.*t44.*4.0;
		t119 = cos(t45);
		t120 = sin(t41);
		t121 = sin(t45);
		t122 = sin(t43);
		t123 = t39.*t40.*t119.*t120.*t121.*t122.*8.0;
		t135 = t42.*2.0;
		t136 = t44.*2.0;
		t124 = t118+t123-t135-t136+1.0;
		t125 = cos(rf2Medy);
		t126 = t30.*t117;
		t127 = t116+t126;
		t128 = t13-t22;
		t129 = t17+t18;
		t130 = t16.*t128;
		t131 = t8.*t12.*t129;
		t132 = t130+t131;
		t133 = t8.*t128;
		t139 = t12.*t16.*t129;
		t134 = t133-t139;
		t138 = t21.*t132;
		t140 = t26.*t134;
		t141 = t138+t140;
		t142 = t21.*t134;
		t144 = t26.*t132;
		t143 = t142-t144;
		t145 = t30.*t141;
		t146 = t145-t35.*t143;
		t148 = t30.*t37;
		t149 = t116+t148;
		t150 = t16.*t21.*t73;
		t151 = t8.*t26.*t73;
		t152 = t150+t151;
		t153 = t8.*t21.*t73;
		t156 = t16.*t26.*t73;
		t154 = t153-t156;
		t155 = t35.*t152;
		t157 = t30.*t154;
		t158 = t155+t157;
		t161 = t24-t31;
		t159 = t26.*t161;
		t160 = t33+t159;
		t170 = t21.*t161;
		t162 = t68-t170;
		t163 = t30.*t162;
		t164 = t163-t35.*t160;
		t165 = t66.*t69;
		t166 = t70.*t73;
		t167 = t165+t166;
		t168 = t66.*t73;
		t169 = t168-t69.*t70;
		out1 = -c.*(q2w3__dt_1_.*(t115.*(t67-t35.*t117)-t70.*t124.*t127+t65.*t66.*t125.*t127)+q2w1__dt_1_.*(t124.*t167-t65.*t125.*t169)+q1flex__dt_1_.*(t91.*(l12z.*t107-l12y.*(t47.*t82.*2.0-t47.*t49.*t82.*4.0-t46.*t80.*t81.*t82.^2.*t83.*4.0+t46.*t51.*t80.*t81.*t83.*4.0))+t103.*(l12y.*t62.*t63-l12z.*t63.*t64)-t74.*t92.*(l12y.*t78+l12z.*t109))+q2w2__dt_1_.*(t69.*t115-t70.*t124.*t149+t65.*t66.*t125.*t149)-q1dev__dt_1_.*(t91.*(t79+t84+l12y.*(t46.*t80.*2.0-t46.*t51.*t80.*4.0-t47.*t80.^2.*t81.*t82.*t83.*4.0+t47.*t49.*t81.*t82.*t83.*4.0))-t74.*t92.*(t95-l12z.*t94+l12y.*t107))-q2w6__dt_1_.*(t124.*(t70.*t146+t9.*t66.*t129)+t115.*(t30.*t143+t35.*t141)-t65.*t125.*(t66.*t146-t9.*t70.*t129))+q2w4__dt_1_.*(t115.*(t30.*t160+t35.*t162)+t70.*t124.*t164-t65.*t66.*t125.*t164)+q2w5__dt_1_.*(t124.*(t15.*t66+t70.*t158)-t115.*(t30.*t152-t35.*t154)+t65.*t125.*(t15.*t70-t66.*t158))-q1sup__dt_1_.*(t91.*(l12y.*(t46.*t47.*t80.*t81.^2.*t82.*4.0-t46.*t47.*t80.*t82.*t83.^2.*4.0)-l12x.*t53.*t61+l12z.*t53.*t62.*t63)+t103.*(l12x.*t63+l12z.*t61.*t62+l12y.*t61.*t64)+t74.*t92.*(-l12x.*t61.*t75+l12z.*t62.*t63.*t75+l12y.*t63.*t64.*t75)))+c.*(l11y-l21y+t91.*(t95+l12y.*(t49.*-2.0-t51.*2.0+t49.*t51.*4.0+t46.*t47.*t80.*t81.*t82.*t83.*8.0+1.0)-l12z.*t94)+t115.*t149-t124.*t169-t103.*(-l12x.*t61+l12z.*t62.*t63+l12y.*t63.*t64)-t65.*t125.*t167+t74.*t92.*(t79+t84-l12y.*t109));
		end
		function out1 = ConsGF_3(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSGF_3
		%    OUT1 = CONSGF_3(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:36
		c = in3(3,:);
		l21 = in3(19,:);
		l22 = in3(20,:);
		l23 = in3(21,:);
		l24 = in3(22,:);
		l25 = in3(23,:);
		l12x = in3(10,:);
		l11z = in3(9,:);
		l12y = in3(11,:);
		l13x = in3(13,:);
		l12z = in3(12,:);
		l13y = in3(14,:);
		l13z = in3(15,:);
		l21z = in3(18,:);
		q1dev__dt_0_ = in2(1,:);
		q1dev__dt_1_ = in2(10,:);
		q1flex__dt_0_ = in2(2,:);
		q1flex__dt_1_ = in2(11,:);
		q1sup__dt_0_ = in2(3,:);
		q1sup__dt_1_ = in2(12,:);
		q2w1__dt_0_ = in2(4,:);
		q2w1__dt_1_ = in2(13,:);
		q2w2__dt_0_ = in2(5,:);
		q2w2__dt_1_ = in2(14,:);
		q2w3__dt_0_ = in2(6,:);
		q2w3__dt_1_ = in2(15,:);
		q2w4__dt_0_ = in2(7,:);
		q2w4__dt_1_ = in2(16,:);
		q2w5__dt_0_ = in2(8,:);
		q2w5__dt_1_ = in2(17,:);
		q2w6__dt_0_ = in2(9,:);
		q2w6__dt_1_ = in2(18,:);
		rf2Medx = in3(33,:);
		rf2Medy = in3(34,:);
		rf1Devx = in3(27,:);
		rf1Devy = in3(28,:);
		riw1 = in3(36,:);
		riw2 = in3(37,:);
		riw3 = in3(38,:);
		riw4 = in3(39,:);
		riw5 = in3(40,:);
		riw6 = in3(41,:);
		t2 = q2w5__dt_0_+riw5;
		t3 = q2w6__dt_0_+riw6;
		t4 = q2w4__dt_0_+riw4;
		t5 = cos(t3);
		t6 = sin(t3);
		t7 = q2w3__dt_0_+riw3;
		t8 = cos(t4);
		t9 = cos(t2);
		t10 = l25+l13z;
		t11 = t9.*t10;
		t12 = sin(t2);
		t13 = l13x.*t5;
		t22 = l13y.*t6;
		t14 = l24+t13-t22;
		t23 = t12.*t14;
		t15 = t11-t23;
		t16 = sin(t4);
		t17 = l13y.*t5;
		t18 = l13x.*t6;
		t19 = l23+t17+t18;
		t20 = q2w2__dt_0_+riw2;
		t21 = cos(t7);
		t24 = t8.*t19;
		t39 = t15.*t16;
		t25 = l22+t24-t39;
		t26 = sin(t7);
		t27 = t8.*t15;
		t28 = t16.*t19;
		t29 = t27+t28;
		t30 = q2w1__dt_0_+riw1;
		t31 = sin(q1dev__dt_0_);
		t32 = sin(q1flex__dt_0_);
		t33 = cos(q1dev__dt_0_);
		t34 = cos(q1flex__dt_0_);
		t35 = sin(q1sup__dt_0_);
		t36 = cos(q1sup__dt_0_);
		t37 = cos(t30);
		t38 = sin(t20);
		t40 = t25.*t26;
		t41 = t21.*t29;
		t42 = t40+t41;
		t43 = t38.*t42;
		t44 = cos(t20);
		t45 = t21.*t25;
		t53 = t26.*t29;
		t46 = l21+t45-t53;
		t67 = t44.*t46;
		t47 = t43-t67;
		t48 = sin(t30);
		t49 = t10.*t12;
		t50 = t9.*t14;
		t51 = t49+t50;
		t52 = cos(rf2Medy);
		t54 = cos(rf1Devy);
		t57 = q1dev__dt_0_./2.0;
		t55 = cos(t57);
		t59 = q1flex__dt_0_./2.0;
		t56 = cos(t59);
		t58 = t55.^2;
		t60 = t56.^2;
		t61 = q1sup__dt_0_./2.0;
		t62 = sin(rf2Medy);
		t63 = t42.*t44;
		t64 = t38.*t46;
		t65 = t63+t64;
		t66 = cos(rf2Medx);
		t68 = sin(rf2Medx);
		t69 = t37.*t47;
		t70 = t69-t48.*t51;
		t71 = t47.*t48;
		t72 = t37.*t51;
		t73 = t71+t72;
		t74 = t16.*t21.*t51;
		t75 = t8.*t26.*t51;
		t76 = t74+t75;
		t77 = t44.*t76;
		t78 = t8.*t21.*t51;
		t82 = t16.*t26.*t51;
		t79 = t78-t82;
		t80 = t38.*t79;
		t81 = t77+t80;
		t83 = sin(rf1Devy);
		t84 = t32.*t33;
		t85 = t84-t31.*t34.*t35;
		t86 = l12x.*t31.*t36;
		t87 = sin(rf1Devx);
		t88 = t31.*t32;
		t89 = t33.*t34.*t35;
		t90 = t88+t89;
		t91 = l12z.*t90;
		t92 = sin(t57);
		t93 = cos(t61);
		t94 = sin(t59);
		t95 = sin(t61);
		t96 = l12x.*t33.*t36;
		t97 = t24-t39;
		t101 = t21.*t97;
		t98 = t53-t101;
		t99 = t26.*(t24-t39);
		t100 = t41+t99;
		t102 = t38.*t98;
		t103 = cos(rf1Devx);
		t104 = t31.*t34;
		t105 = t104-t32.*t33.*t35;
		t106 = t33.*t34;
		t107 = t31.*t32.*t35;
		t108 = t106+t107;
		t109 = t45-t53;
		t110 = t38.*t109;
		t111 = t63+t110;
		t112 = t13-t22;
		t113 = t17+t18;
		t114 = t16.*t112;
		t115 = t8.*t12.*t113;
		t116 = t114+t115;
		t117 = t8.*t112;
		t119 = t12.*t16.*t113;
		t118 = t117-t119;
		t120 = t21.*t118;
		t127 = t26.*t116;
		t121 = t120-t127;
		t122 = t44.*t121;
		t123 = t21.*t116;
		t124 = t26.*t118;
		t125 = t123+t124;
		t126 = t122-t38.*t125;
		out1 = c.*(l11z-l21z+t62.*t73-t83.*(t91+t96-l12y.*t105)+t54.*t103.*(-l12x.*t35+l12y.*t32.*t36+l12z.*t34.*t36)+t54.*t87.*(t86+l12y.*(t58.*-2.0-t60.*2.0+t58.*t60.*4.0+t55.*t56.*t92.*t93.*t94.*t95.*8.0+1.0)-l12z.*t85)-t52.*t65.*t66+t52.*t68.*t70)-c.*(-q2w1__dt_1_.*(t62.*t70-t52.*t68.*t73)+q2w6__dt_1_.*(t62.*(t48.*t126+t9.*t37.*t113)+t52.*t66.*(t38.*t121+t44.*t125)+t52.*t68.*(t37.*t126-t9.*t48.*t113))+q2w4__dt_1_.*(t48.*t62.*(t102-t44.*(t41+t26.*t97))-t52.*t66.*(t38.*t100+t44.*t98)+t37.*t52.*t68.*(t102-t44.*t100))+q2w5__dt_1_.*(-t62.*(t15.*t37-t48.*t81)+t52.*t68.*(t15.*t48+t37.*t81)+t52.*t66.*(t38.*t76-t44.*t79))+q1sup__dt_1_.*(t83.*(-l12x.*t33.*t35+l12y.*t32.*t33.*t36+l12z.*t33.*t34.*t36)+t54.*t103.*(l12x.*t36+l12y.*t32.*t35+l12z.*t34.*t35)-t54.*t87.*(l12y.*(t55.*t56.*t92.*t93.^2.*t94.*4.0-t55.*t56.*t92.*t94.*t95.^2.*4.0)-l12x.*t31.*t35+l12z.*t31.*t34.*t36))-q2w3__dt_1_.*(t52.*t66.*(t43-t44.*t109)+t48.*t62.*t111+t37.*t52.*t68.*t111)+q1flex__dt_1_.*(t83.*(l12y.*t90+l12z.*t105)+t54.*t103.*(l12z.*t32.*t36-l12y.*t34.*t36)+t54.*t87.*(l12z.*t108-l12y.*(t56.*t94.*2.0-t56.*t58.*t94.*4.0-t55.*t92.*t93.*t94.^2.*t95.*4.0+t55.*t60.*t92.*t93.*t95.*4.0)))-q1dev__dt_1_.*(t83.*(t86-l12z.*t85+l12y.*t108)+t54.*t87.*(t91+t96+l12y.*(t55.*t92.*2.0-t55.*t60.*t92.*4.0-t56.*t92.^2.*t93.*t94.*t95.*4.0+t56.*t58.*t93.*t94.*t95.*4.0)))-q2w2__dt_1_.*(t47.*t52.*t66+t48.*t62.*t65+t37.*t52.*t65.*t68));
		end
		function out1 = ConsGF_4(t,in2,in3,in4,in5,in6)
		%CONSGF_4
		%    OUT1 = CONSGF_4(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:39
		ax__FRAME_5__dt_1_ = in6(7);
		ax__FRAME_20__dt_1_ = in6(20);
		ay__FRAME_5__dt_1_ = in6(8);
		ay__FRAME_20__dt_1_ = in6(21);
		az__FRAME_5__dt_1_ = in6(9);
		az__FRAME_20__dt_1_ = in6(22);
		c = in3(3,:);
		rqw__FRAME_5__dt_0_ = in6(10);
		rqw__FRAME_20__dt_0_ = in6(23);
		rqx__FRAME_5__dt_0_ = in6(11);
		rqx__FRAME_20__dt_0_ = in6(24);
		rqy__FRAME_5__dt_0_ = in6(12);
		rqy__FRAME_20__dt_0_ = in6(25);
		rqz__FRAME_5__dt_0_ = in6(13);
		rqz__FRAME_20__dt_0_ = in6(26);
		out1 = -c.*(rqw__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_-rqw__FRAME_20__dt_0_.*rqx__FRAME_5__dt_0_+rqy__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_-rqy__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)+c.*(rqw__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0)+rqy__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0-(ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0)+rqx__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0)-rqz__FRAME_20__dt_0_.*(ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_.*(-1.0./2.0)+(ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0)-rqw__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0)-rqy__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0-(ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0)-rqx__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0)+rqz__FRAME_5__dt_0_.*(ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_.*(-1.0./2.0)+(ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0));
		end
		function out1 = ConsGF_5(t,in2,in3,in4,in5,in6)
		%CONSGF_5
		%    OUT1 = CONSGF_5(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:42
		ax__FRAME_5__dt_1_ = in6(7);
		ax__FRAME_20__dt_1_ = in6(20);
		ay__FRAME_5__dt_1_ = in6(8);
		ay__FRAME_20__dt_1_ = in6(21);
		az__FRAME_5__dt_1_ = in6(9);
		az__FRAME_20__dt_1_ = in6(22);
		c = in3(3,:);
		rqw__FRAME_5__dt_0_ = in6(10);
		rqw__FRAME_20__dt_0_ = in6(23);
		rqx__FRAME_5__dt_0_ = in6(11);
		rqx__FRAME_20__dt_0_ = in6(24);
		rqy__FRAME_5__dt_0_ = in6(12);
		rqy__FRAME_20__dt_0_ = in6(25);
		rqz__FRAME_5__dt_0_ = in6(13);
		rqz__FRAME_20__dt_0_ = in6(26);
		out1 = -c.*(rqw__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_-rqw__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_-rqx__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_+rqx__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)+c.*(rqw__FRAME_20__dt_0_.*(ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_.*(-1.0./2.0)+(ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0)-rqx__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0-(ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0)-rqw__FRAME_5__dt_0_.*(ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_.*(-1.0./2.0)+(ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0)+rqx__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0-(ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0)+rqy__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0)+rqz__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0)-rqy__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0)-rqz__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0));
		end
		function out1 = ConsGF_6(t,in2,in3,in4,in5,in6)
		%CONSGF_6
		%    OUT1 = CONSGF_6(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:45
		ax__FRAME_5__dt_1_ = in6(7);
		ax__FRAME_20__dt_1_ = in6(20);
		ay__FRAME_5__dt_1_ = in6(8);
		ay__FRAME_20__dt_1_ = in6(21);
		az__FRAME_5__dt_1_ = in6(9);
		az__FRAME_20__dt_1_ = in6(22);
		c = in3(3,:);
		rqw__FRAME_5__dt_0_ = in6(10);
		rqw__FRAME_20__dt_0_ = in6(23);
		rqx__FRAME_5__dt_0_ = in6(11);
		rqx__FRAME_20__dt_0_ = in6(24);
		rqy__FRAME_5__dt_0_ = in6(12);
		rqy__FRAME_20__dt_0_ = in6(25);
		rqz__FRAME_5__dt_0_ = in6(13);
		rqz__FRAME_20__dt_0_ = in6(26);
		out1 = -c.*(rqw__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_-rqw__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_+rqx__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_-rqx__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_)+c.*(rqw__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0-(ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0)+rqx__FRAME_20__dt_0_.*(ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_.*(-1.0./2.0)+(ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0)-rqw__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0-(ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0)-rqy__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0)-rqx__FRAME_5__dt_0_.*(ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_.*(-1.0./2.0)+(ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0)+rqz__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0)+rqy__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0)-rqz__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0));
		end
		function out1 = ConsGF_7(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSGF_7
		%    OUT1 = CONSGF_7(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:20:20
		c = in3(3,:);
		q1dev__dt_0_ = in2(1,:);
		q1dev__dt_1_ = in2(10,:);
		q1flex__dt_0_ = in2(2,:);
		q1flex__dt_1_ = in2(11,:);
		q1sup__dt_0_ = in2(3,:);
		q1sup__dt_1_ = in2(12,:);
		xx3__dt_0_ = in5(1,:);
		zz3__dt_0_ = in5(3,:);
		t2 = q1flex__dt_0_./2.0;
		t3 = q1dev__dt_0_./2.0;
		t4 = q1sup__dt_0_./2.0;
		t5 = cos(t2);
		t6 = cos(t4);
		t7 = sin(t3);
		t8 = cos(t3);
		t9 = sin(t2);
		t10 = sin(t4);
		t11 = t6.^2;
		t12 = t8.^2;
		t13 = t5.^2;
		t14 = t11.^2;
		t15 = cos(q1sup__dt_0_);
		t16 = sin(q1sup__dt_0_);
		t17 = xx3__dt_0_./2.0;
		t18 = zz3__dt_0_./2.0;
		t19 = t11.*t12.*2.0;
		t20 = t11.*t13.*2.0;
		t21 = t7.*t9.*t10.*2.0;
		t22 = t5.*t6.*t8.*t11.*2.0;
		t23 = t11.*-3.0-t12-t13+t14+t19+t20+t21+t22-t7.*t9.*t10.*t11.*2.0+2.0;
		t24 = 1.0./t23;
		t25 = sin(q1dev__dt_0_);
		t26 = cos(q1dev__dt_0_);
		t27 = cos(q1flex__dt_0_);
		t28 = sin(q1flex__dt_0_);
		t29 = sin(t17);
		t30 = sin(t18);
		out1 = -c.*(-q1dev__dt_1_.*((t5.*t7.*t10)./2.0-(t6.*t8.*t9)./2.0)+q1flex__dt_1_.*((t5.*t6.*t7)./2.0-(t8.*t9.*t10)./2.0)+q1sup__dt_1_.*((t5.*t6.*t8)./2.0-(t7.*t9.*t10)./2.0)+(t24.*t29.*cos(t18).*(q1dev__dt_1_.*t14.*-2.0+q1dev__dt_1_.*t15+(q1flex__dt_1_.*t16)./2.0+(q1sup__dt_1_.*t28)./2.0-(q1dev__dt_1_.*t15.*t26)./2.0-(q1dev__dt_1_.*t15.*t27)./2.0+(q1sup__dt_1_.*t16.*t25)./2.0+q1dev__dt_1_.*t5.*t6.*t8-(q1dev__dt_1_.*t7.*t9.*t10)./2.0+q1flex__dt_1_.*t6.*t7.*t9+(q1flex__dt_1_.*t5.*t8.*t10)./2.0+q1sup__dt_1_.*t5.*t7.*t10.*(3.0./2.0)+q1sup__dt_1_.*t6.*t8.*t9.*2.0-q1dev__dt_1_.*t5.*t6.*t8.*t11.*3.0+q1dev__dt_1_.*t7.*t9.*t10.*t15.*(3.0./2.0)-q1flex__dt_1_.*t6.*t7.*t9.*t11+(q1flex__dt_1_.*t5.*t8.*t10.*t15)./2.0-q1sup__dt_1_.*t6.*t8.*t9.*t11+(q1sup__dt_1_.*t5.*t7.*t10.*t15)./2.0))./4.0+(t24.*t30.*cos(t17).*((q1dev__dt_1_.*t16)./2.0-q1flex__dt_1_.*t14.*2.0+q1flex__dt_1_.*t15+(q1sup__dt_1_.*t25)./2.0-(q1flex__dt_1_.*t15.*t26)./2.0-(q1flex__dt_1_.*t15.*t27)./2.0+(q1sup__dt_1_.*t16.*t28)./2.0+q1dev__dt_1_.*t6.*t7.*t9+(q1dev__dt_1_.*t5.*t8.*t10)./2.0+q1flex__dt_1_.*t5.*t6.*t8-(q1flex__dt_1_.*t7.*t9.*t10)./2.0+q1sup__dt_1_.*t5.*t6.*t7.*2.0+q1sup__dt_1_.*t8.*t9.*t10.*(3.0./2.0)-q1dev__dt_1_.*t6.*t7.*t9.*t11+(q1dev__dt_1_.*t5.*t8.*t10.*t15)./2.0-q1flex__dt_1_.*t5.*t6.*t8.*t11.*3.0+q1flex__dt_1_.*t7.*t9.*t10.*t15.*(3.0./2.0)-q1sup__dt_1_.*t5.*t6.*t7.*t11+(q1sup__dt_1_.*t8.*t9.*t10.*t15)./2.0))./4.0)-c.*(t29.*t30.*(-1.0./2.0)+t6.*t7.*t9+t5.*t8.*t10);
		end
		function out1 = ConsFVJac_1(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSFVJAC_1
		%    OUT1 = CONSFVJAC_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:17:58
		out1 = 0.0;
		end
		function out1 = ConsFVJac_2(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSFVJAC_2
		%    OUT1 = CONSFVJAC_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:18:55
		out1 = 0.0;
		end
		function out1 = ConsFVJac_3(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSFVJAC_3
		%    OUT1 = CONSFVJAC_3(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:37
		out1 = 0.0;
		end
		function out1 = ConsFVJac_4(t,in2,in3,in4,in5,in6)
		%CONSFVJAC_4
		%    OUT1 = CONSFVJAC_4(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:40
		rqw__FRAME_5__dt_0_ = in6(10);
		rqw__FRAME_20__dt_0_ = in6(23);
		rqx__FRAME_5__dt_0_ = in6(11);
		rqx__FRAME_20__dt_0_ = in6(24);
		rqy__FRAME_5__dt_0_ = in6(12);
		rqy__FRAME_20__dt_0_ = in6(25);
		rqz__FRAME_5__dt_0_ = in6(13);
		rqz__FRAME_20__dt_0_ = in6(26);
		t2 = (rqw__FRAME_5__dt_0_.*rqw__FRAME_20__dt_0_)./2.0;
		t3 = (rqx__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_)./2.0;
		t4 = (rqy__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
		t5 = (rqz__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
		t6 = (rqw__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0;
		t7 = (rqx__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_)./2.0;
		t8 = t6+t7-(rqw__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0-(rqx__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
		t9 = (rqw__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
		t10 = (rqx__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0;
		t11 = t9+t10-(rqw__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_)./2.0-(rqx__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,t2+t3+t4+t5,-t2-t3-t4-t5,t8,t8,t11,t11],[2,6]);
		end
		function out1 = ConsFVJac_5(t,in2,in3,in4,in5,in6)
		%CONSFVJAC_5
		%    OUT1 = CONSFVJAC_5(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:43
		rqw__FRAME_5__dt_0_ = in6(10);
		rqw__FRAME_20__dt_0_ = in6(23);
		rqx__FRAME_5__dt_0_ = in6(11);
		rqx__FRAME_20__dt_0_ = in6(24);
		rqy__FRAME_5__dt_0_ = in6(12);
		rqy__FRAME_20__dt_0_ = in6(25);
		rqz__FRAME_5__dt_0_ = in6(13);
		rqz__FRAME_20__dt_0_ = in6(26);
		t2 = (rqw__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
		t3 = (rqx__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
		t4 = t2+t3-(rqw__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0-(rqx__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_)./2.0;
		t5 = (rqw__FRAME_5__dt_0_.*rqw__FRAME_20__dt_0_)./2.0;
		t6 = (rqx__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_)./2.0;
		t7 = (rqy__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
		t8 = (rqz__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
		t9 = (rqw__FRAME_20__dt_0_.*rqx__FRAME_5__dt_0_)./2.0;
		t10 = (rqy__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0;
		t11 = t9+t10-(rqw__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_)./2.0-(rqy__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,t4,t4,t5+t6+t7+t8,-t5-t6-t7-t8,t11,t11],[2,6]);
		end
		function out1 = ConsFVJac_6(t,in2,in3,in4,in5,in6)
		%CONSFVJAC_6
		%    OUT1 = CONSFVJAC_6(T,IN2,IN3,IN4,IN5,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:19:46
		rqw__FRAME_5__dt_0_ = in6(10);
		rqw__FRAME_20__dt_0_ = in6(23);
		rqx__FRAME_5__dt_0_ = in6(11);
		rqx__FRAME_20__dt_0_ = in6(24);
		rqy__FRAME_5__dt_0_ = in6(12);
		rqy__FRAME_20__dt_0_ = in6(25);
		rqz__FRAME_5__dt_0_ = in6(13);
		rqz__FRAME_20__dt_0_ = in6(26);
		t2 = (rqw__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_)./2.0;
		t3 = (rqx__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
		t4 = t2+t3-(rqw__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0-(rqx__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0;
		t5 = (rqw__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_)./2.0;
		t6 = (rqy__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
		t7 = t5+t6-(rqw__FRAME_20__dt_0_.*rqx__FRAME_5__dt_0_)./2.0-(rqy__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0;
		t8 = (rqw__FRAME_5__dt_0_.*rqw__FRAME_20__dt_0_)./2.0;
		t9 = (rqx__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_)./2.0;
		t10 = (rqy__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
		t11 = (rqz__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,t4,t4,t7,t7,t8+t9+t10+t11,-t8-t9-t10-t11],[2,6]);
		end
		function out1 = ConsFVJac_7(t,in2,in3,in4,in5,SUBSVECTOR__)
		%CONSFVJAC_7
		%    OUT1 = CONSFVJAC_7(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    04-Jan-2021 17:20:21
		out1 = 0.0;
		end
	end

end

