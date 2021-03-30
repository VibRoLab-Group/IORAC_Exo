function [TF,Vel,Cor,Jac,TransDis,Quat,MM,GFI,GFF,GFU,JacU,JacCons,CorCons,GFCons,CenMat,CorMatLeft,CorMatRight]=System_Arm4D(t,q,p,u,s)
    
    [TF,Vel,Cor,Jac,TransDis,Quat]=numKinematics_Arm4D(t,q,p,u,s);
    [MM,GFI,CenMat,CorMatLeft,CorMatRight]=Inertial_Arm4D(t,q,p,u,s,TF,Vel,Cor,Jac);
    [GFF]=Force_Arm4D(t,q,p,u,s,TF,TransDis,Vel,Jac,Quat);
    [GFU,JacU]=Input_Arm4D(t,q,p,u,s,TF,TransDis,Vel,Jac,Quat);
    [JacCons,CorCons,GFCons]=Constraint_Arm4D(t,q,p,u,s,TransDis,Vel,Cor,Jac,Quat);

	function [frameTF,frameVel,frameCorAcc,frameJacobian,frameTransDis,frameRotQuat]=numKinematics_Arm4D(t,q,p,u,s)
	    LinkTF=linkTF(t,q,p,u,s);
	    LinkVel=linkVel(t,q,p,u,s);
	    LinkCorAcc=linkCorAcc(t,q,p,u,s);
	    LinkJacobian=linkJacobian(t,q,p,u,s);
	    LinkRotQuat=linkRotQuat(t,q,p,u,s);
	    FrameNum=19;
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
			        framePath=[1  2  3  8];
			    case 9
			        framePath=[1  2  3  4  9];
			    case 10
			        framePath=[1   2   3   4   5  10];
			    case 11
			        framePath=[1   2   3   4   5   6  11];
			    case 12
			        framePath=[1   2   3   4   5  10  12];
			    case 13
			        framePath=[1   2   3   4   5  10  13];
			    case 14
			        framePath=[1   2   3   4   5  10  14];
			    case 15
			        framePath=[1   2   3   4   5  10  15];
			    case 16
			        framePath=[1   2   3   4   5   6  11  16];
			    case 17
			        framePath=[1   2   3   4   5   6  11  17];
			    case 18
			        framePath=[1   2   3   4   5   6  11  18];
			    case 19
			        framePath=[1   2   3   4   5   6  11  19];
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
		function out1 = linkTF(t,in2,in3,in4,NHSIGNAL)
		%LINKTF
		%    OUT1 = LINKTF(T,IN2,IN3,IN4,NHSIGNAL)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:57:01
		ang0x = in3(56,:);
		ang0y = in3(57,:);
		ang0z = in3(58,:);
		link0x = in3(29,:);
		link0y = in3(30,:);
		link1x = in3(32,:);
		link0z = in3(31,:);
		link1y = in3(33,:);
		link2x = in3(35,:);
		link1z = in3(34,:);
		link2y = in3(36,:);
		link3x = in3(38,:);
		link2z = in3(37,:);
		link3y = in3(39,:);
		link4x = in3(41,:);
		link3z = in3(40,:);
		link4y = in3(42,:);
		link4z = in3(43,:);
		linkCOM1x = in3(44,:);
		linkCOM1y = in3(45,:);
		linkCOM2x = in3(47,:);
		linkCOM1z = in3(46,:);
		linkCOM2y = in3(48,:);
		linkCOM3x = in3(50,:);
		linkCOM2z = in3(49,:);
		linkCOM3y = in3(51,:);
		linkCOM4x = in3(53,:);
		linkCOM3z = in3(52,:);
		linkCOM4y = in3(54,:);
		linkCOM4z = in3(55,:);
		qJ1__dt_0_ = in2(1,:);
		qJ2__dt_0_ = in2(2,:);
		qJ3__dt_0_ = in2(3,:);
		qJ4__dt_0_ = in2(4,:);
		rTetra = in3(69,:);
		t2 = (cos(ang0z));
		t3 = (sin(ang0x));
		t4 = (sin(ang0z));
		t5 = (cos(ang0x));
		t6 = (sin(ang0y));
		t7 = (sqrt(3.0));
		t8 = (cos(ang0y));
		t9 = (cos(qJ1__dt_0_));
		t10 = (sin(qJ1__dt_0_));
		t11 = (sqrt(6.0));
		t12 = ((rTetra.*t11)./2.0);
		t13 = (sin(qJ2__dt_0_));
		t14 = (cos(qJ2__dt_0_));
		t15 = (sin(qJ3__dt_0_));
		t16 = (cos(qJ3__dt_0_));
		t17 = (sin(qJ4__dt_0_));
		t18 = (cos(qJ4__dt_0_));
		t19 = (rTetra.*t7.*(2.0./3.0));
		outSize=[4,76];
		elementRow=[1,2,3,4,1,2,3,1,2,3,1,2,3,4,1,2,3,1,2,1,2,3,1,2,3,4,1,2,3,2,3,1,2,3,4,1,3,2,1,3,1,2,3,4,1,2,3,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,2,3,4,1,2,3,2,4,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,2,3,4,1,2,3,2,4];
		elementCol=[1,2,3,4,5,5,5,6,6,6,7,7,7,8,9,9,9,10,10,11,11,11,12,12,12,12,13,14,14,15,15,16,16,16,16,17,17,18,19,19,20,20,20,20,21,22,22,23,23,24,24,24,24,25,26,27,28,28,28,28,29,30,31,32,32,32,32,33,34,35,36,36,36,36,37,38,39,40,40,40,40,41,42,43,44,44,44,44,45,46,47,48,48,48,48,49,50,51,52,52,52,52,53,54,55,56,56,56,57,58,59,60,60,61,62,63,64,64,64,64,65,66,67,68,68,68,68,69,70,71,72,72,72,73,74,75,76,76];
		elementList=[1.0,1.0,1.0,1.0,t2.*t8,t4.*t8,-t6,-t4.*t5+t2.*t3.*t6,t2.*t5+t3.*t4.*t6,t3.*t8,t3.*t4+t2.*t5.*t6,-t2.*t3+t4.*t5.*t6,t5.*t8,1.0,(t7.*t9)./2.0,t9.*(-1.0./2.0),-t10,1.0./2.0,t7./2.0,(t7.*t10)./2.0,t10.*(-1.0./2.0),t9,link0x,link0y,link0z,1.0,1.0,t14,t13,-t13,t14,link1x,link1y,link1z,1.0,t16,-t15,1.0,t15,t16,link2x,link2y,link2z,1.0,1.0,t18,t17,-t17,t18,link3x,link3y,link3z,1.0,1.0,1.0,1.0,link4x,link4y,link4z,1.0,1.0,1.0,1.0,linkCOM1x,linkCOM1y,linkCOM1z,1.0,1.0,1.0,1.0,linkCOM2x,linkCOM2y,linkCOM2z,1.0,1.0,1.0,1.0,linkCOM3x,linkCOM3y,linkCOM3z,1.0,1.0,1.0,1.0,linkCOM4x,linkCOM4y,linkCOM4z,1.0,1.0,1.0,1.0,-rTetra,rTetra.*t11.*(-1.0./6.0),rTetra.*t7.*(-1.0./3.0),1.0,1.0,1.0,1.0,rTetra,rTetra.*t11.*(-1.0./6.0),rTetra.*t7.*(-1.0./3.0),1.0,1.0,1.0,1.0,rTetra.*t11.*(-1.0./6.0),t19,1.0,1.0,1.0,1.0,t12,1.0,1.0,1.0,1.0,-rTetra,rTetra.*t11.*(-1.0./6.0),rTetra.*t7.*(-1.0./3.0),1.0,1.0,1.0,1.0,rTetra,rTetra.*t11.*(-1.0./6.0),rTetra.*t7.*(-1.0./3.0),1.0,1.0,1.0,1.0,rTetra.*t11.*(-1.0./6.0),t19,1.0,1.0,1.0,1.0,t12,1.0];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = linkVel(t,in2,in3,in4,NHSIGNAL)
		%LINKVEL
		%    OUT1 = LINKVEL(T,IN2,IN3,IN4,NHSIGNAL)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:57:00
		qJ1__dt_1_ = in2(5,:);
		qJ2__dt_1_ = in2(6,:);
		qJ3__dt_1_ = in2(7,:);
		qJ4__dt_1_ = in2(8,:);
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,qJ1__dt_1_./2.0,(sqrt(3.0).*qJ1__dt_1_)./2.0,0.0,0.0,0.0,0.0,qJ2__dt_1_,0.0,0.0,0.0,0.0,0.0,0.0,qJ3__dt_1_,0.0,0.0,0.0,0.0,qJ4__dt_1_,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,19]);
		end
		function out1 = linkRotQuat(t,in2,in3,in4,NHSIGNAL)
		%LINKROTQUAT
		%    OUT1 = LINKROTQUAT(T,IN2,IN3,IN4,NHSIGNAL)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:57:00
		ang0x = in3(56,:);
		ang0y = in3(57,:);
		ang0z = in3(58,:);
		qJ1__dt_0_ = in2(1,:);
		qJ2__dt_0_ = in2(2,:);
		qJ3__dt_0_ = in2(3,:);
		qJ4__dt_0_ = in2(4,:);
		t2 = ang0x./2.0;
		t3 = ang0y./2.0;
		t4 = ang0z./2.0;
		t5 = cos(t3);
		t6 = cos(t4);
		t7 = sin(t2);
		t8 = cos(t2);
		t9 = sin(t3);
		t10 = sin(t4);
		t11 = qJ1__dt_0_./2.0;
		t12 = sqrt(2.0);
		t13 = sqrt(6.0);
		t14 = qJ2__dt_0_./2.0;
		t15 = qJ4__dt_0_./2.0;
		t16 = sin(t11);
		t17 = t12+t13;
		t18 = qJ3__dt_0_./2.0;
		t19 = cos(t11);
		t20 = t12-t13;
		out1 = reshape([1.0,0.0,0.0,0.0,t5.*t6.*t8+t7.*t9.*t10,t5.*t6.*t7-t8.*t9.*t10,t5.*t7.*t10+t6.*t8.*t9,-t6.*t7.*t9+t5.*t8.*t10,(t17.*t19)./4.0,t16.*t20.*(-1.0./4.0),(t16.*t17)./4.0,(t19.*t20)./4.0,cos(t14),sin(t14),0.0,0.0,cos(t18),0.0,sin(t18),0.0,cos(t15),sin(t15),0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0],[4,19]);
		end
		function out1 = linkCorAcc(t,in2,in3,in4,NHSIGNAL)
		%LINKCORACC
		%    OUT1 = LINKCORACC(T,IN2,IN3,IN4,NHSIGNAL)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:57:01
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,19]);
		end
		function out1 = linkJacobian(t,in2,in3,in4,NHSIGNAL)
		%LINKJACOBIAN
		%    OUT1 = LINKJACOBIAN(T,IN2,IN3,IN4,NHSIGNAL)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:57:02
		outSize=[6,76];
		elementRow=[4,5,4,5,4];
		elementCol=[9,9,14,19,24];
		elementList=[1.0./2.0,sqrt(3.0)./2.0,1.0,1.0,1.0];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
	end

	function [inertM,inertGF,inertCenMat,inertCorMatLeft,inertCorMatRight]=Inertial_Arm4D(t,q,p,u,s,TF_Global,Vel_Global,Cor_Global,Jac_Global)
	    BaseFrameList=[3  4  5  6  5  6];
	    MassList=Mass_Arm4D(t,q,p,s,u);
	    MomentList=reshape(Moment_Arm4D(t,q,p,s,u),[3 3 numel(MassList)]);
	    
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
	        
	        % thisInertCorMatLeft=[m*vJacobian;I*wJacobian].';
	        % thisInertCorMatRight=[vJacobian;wJacobian]; %Take the time derivative of this to calculate coriolis numerically
	        inertCorMatLeft(:,6*bodyNum+(-5:0))=[m*vJacobian;I*wJacobian].';
	        inertCorMatRight(6*bodyNum+(-5:0),:)=[vJacobian;wJacobian];
	        inertCenMat(:,:,bodyNum)=- wJacobian.'* skew3_InertDynamicsOnly_(I*w) * wJacobian; % On the same side of Mqddot
	        % inertM(:,:,bodyNum)=thisInertCorMatLeft*thisInertCorMatRight;
	        % inertGF(:,bodyNum)= - thisInertCorMatLeft*[a;alpha] - wJacobian.'* cross(w,I*w); % on the different side of Mqddot
	        inertM(:,:,bodyNum)=vJacobian.'*m*vJacobian+wJacobian.'*I*wJacobian;
	        inertGF(:,bodyNum)=-vJacobian.'*m*a-wJacobian.'*(I*alpha+cross(w,I*w)); % on the different side of Mqddot
	    end
	    function output_=skew3_InertDynamicsOnly_(w_)
	        output_=[0 -w_(3) w_(2) ; w_(3) 0 -w_(1) ; -w_(2) w_(1) 0 ];
	    end
		function out1 = Mass_Arm4D(t,in2,in3,in4,NHSIGNAL)
		%MASS_ARM4D
		%    OUT1 = MASS_ARM4D(T,IN2,IN3,IN4,NHSIGNAL)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:16
		link1m = in3(1,:);
		link2m = in3(8,:);
		link3m = in3(15,:);
		link4m = in3(22,:);
		linkUkn3m = in3(59,:);
		linkUkn4m = in3(63,:);
		out1 = [link1m,link2m,link3m,link4m,linkUkn3m,linkUkn4m];
		end
		function out1 = Moment_Arm4D(t,in2,in3,in4,NHSIGNAL)
		%MOMENT_ARM4D
		%    OUT1 = MOMENT_ARM4D(T,IN2,IN3,IN4,NHSIGNAL)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:16
		link1i11 = in3(2,:);
		link1i12 = in3(3,:);
		link1i13 = in3(4,:);
		link1i22 = in3(5,:);
		link1i23 = in3(6,:);
		link1i33 = in3(7,:);
		link2i11 = in3(9,:);
		link2i12 = in3(10,:);
		link2i13 = in3(11,:);
		link2i22 = in3(12,:);
		link2i23 = in3(13,:);
		link2i33 = in3(14,:);
		link3i11 = in3(16,:);
		link3i12 = in3(17,:);
		link3i13 = in3(18,:);
		link3i22 = in3(19,:);
		link3i23 = in3(20,:);
		link3i33 = in3(21,:);
		link4i11 = in3(23,:);
		link4i12 = in3(24,:);
		link4i13 = in3(25,:);
		link4i22 = in3(26,:);
		link4i23 = in3(27,:);
		link4i33 = in3(28,:);
		linkUkn3i11 = in3(60,:);
		linkUkn3i22 = in3(61,:);
		linkUkn3i33 = in3(62,:);
		linkUkn4i11 = in3(64,:);
		linkUkn4i22 = in3(65,:);
		linkUkn4i33 = in3(66,:);
		outSize=[3,18];
		elementRow=[1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3];
		elementCol=[1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,14,15,16,17,18];
		elementList=[link1i11,link1i12,link1i13,link1i12,link1i22,link1i23,link1i13,link1i23,link1i33,link2i11,link2i12,link2i13,link2i12,link2i22,link2i23,link2i13,link2i23,link2i33,link3i11,link3i12,link3i13,link3i12,link3i22,link3i23,link3i13,link3i23,link3i33,link4i11,link4i12,link4i13,link4i12,link4i22,link4i23,link4i13,link4i23,link4i33,linkUkn3i11,linkUkn3i22,linkUkn3i33,linkUkn4i11,linkUkn4i22,linkUkn4i33];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
	end

	function [CollectGF]=Force_Arm4D(t,q,p,u,s,TF_Global,TransDis_Global,Vel_Global,Jac_Global,Quat_Global)
	    CollectGF=zeros(numel(q)/2,1);
	    ForceNum=5;
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
			        SubFrame=[8];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[8];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_2(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			    case 3
			        SubFrame=[9];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[9];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_3(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			    case 4
			        SubFrame=[10];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[10];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_4(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			    case 5
			        SubFrame=[11];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[11];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_5(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			end
	        CollectGF(:,FCount)=SubGF;
	    end
		function out1 = ForceJac_1(t,in2,in3,in4,NHSIGNAL,SUBSVECTOR__)
		%FORCEJAC_1
		%    OUT1 = FORCEJAC_1(T,IN2,IN3,IN4,NHSIGNAL,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:17
		outSize=[4,4];
		elementRow=[1,2,3,4];
		elementCol=[1,2,3,4];
		elementList=[1.0,1.0,1.0,1.0];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = ForceJac_2(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_2
		%    OUT1 = FORCEJAC_2(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:18
		out1 = 0.0;
		end
		function out1 = ForceJac_3(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_3
		%    OUT1 = FORCEJAC_3(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:19
		out1 = 0.0;
		end
		function out1 = ForceJac_4(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_4
		%    OUT1 = FORCEJAC_4(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:20
		out1 = 0.0;
		end
		function out1 = ForceJac_5(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_5
		%    OUT1 = FORCEJAC_5(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:21
		out1 = 0.0;
		end
		function out1 = ForceEff_1(t,in2,in3,in4,NHSIGNAL,SUBSVECTOR__)
		%FORCEEFF_1
		%    OUT1 = FORCEEFF_1(T,IN2,IN3,IN4,NHSIGNAL,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:17
		cJoint = in3(68,:);
		qJ1__dt_1_ = in2(5,:);
		qJ2__dt_1_ = in2(6,:);
		qJ3__dt_1_ = in2(7,:);
		qJ4__dt_1_ = in2(8,:);
		out1 = [-cJoint.*qJ1__dt_1_;-cJoint.*qJ2__dt_1_;-cJoint.*qJ3__dt_1_;-cJoint.*qJ4__dt_1_];
		end
		function out1 = ForceEff_2(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEEFF_2
		%    OUT1 = FORCEEFF_2(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:18
		gAcc = in3(67,:);
		link1m = in3(1,:);
		out1 = [0.0;0.0;-gAcc.*link1m];
		end
		function out1 = ForceEff_3(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEEFF_3
		%    OUT1 = FORCEEFF_3(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:19
		gAcc = in3(67,:);
		link2m = in3(8,:);
		out1 = [0.0;0.0;-gAcc.*link2m];
		end
		function out1 = ForceEff_4(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEEFF_4
		%    OUT1 = FORCEEFF_4(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:20
		gAcc = in3(67,:);
		link3m = in3(15,:);
		out1 = [0.0;0.0;-gAcc.*link3m];
		end
		function out1 = ForceEff_5(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEEFF_5
		%    OUT1 = FORCEEFF_5(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:21
		gAcc = in3(67,:);
		link4m = in3(22,:);
		out1 = [0.0;0.0;-gAcc.*link4m];
		end
		function out1 = ForceFVJac_1(t,in2,in3,in4,NHSIGNAL,SUBSVECTOR__)
		%FORCEFVJAC_1
		%    OUT1 = FORCEFVJAC_1(T,IN2,IN3,IN4,NHSIGNAL,SUBSVECTOR__)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:17
		out1 = 0.0;
		end
		function out1 = ForceFVJac_2(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEFVJAC_2
		%    OUT1 = FORCEFVJAC_2(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:18
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = ForceFVJac_3(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEFVJAC_3
		%    OUT1 = FORCEFVJAC_3(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:19
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = ForceFVJac_4(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEFVJAC_4
		%    OUT1 = FORCEFVJAC_4(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:20
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = ForceFVJac_5(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEFVJAC_5
		%    OUT1 = FORCEFVJAC_5(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:22
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
	end

	function [CollectGF,CollectInputJac]=Input_Arm4D(t,q,p,u,s,TF_Global,TransDis_Global,Vel_Global,Jac_Global,Quat_Global)
	    CollectGF=zeros(numel(q)/2,1);
	    CollectInputJac=zeros(numel(q)/2,numel(u));
	    ForceNum=15;
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
			        SubFrame=[12];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[12];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_1(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_1(t,q,p,u,s,SubSubs);
			    case 2
			        SubFrame=[13];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[13];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_2(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_2(t,q,p,u,s,SubSubs);
			    case 3
			        SubFrame=[14];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[14];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_3(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_3(t,q,p,u,s,SubSubs);
			    case 4
			        SubFrame=[15];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[15];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_4(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_4(t,q,p,u,s,SubSubs);
			    case 5
			        SubFrame=[16];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[16];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_5(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_5(t,q,p,u,s,SubSubs);
			    case 6
			        SubFrame=[17];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[17];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_6(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_6(t,q,p,u,s,SubSubs);
			    case 7
			        SubFrame=[18];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[18];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_7(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_7(t,q,p,u,s,SubSubs);
			    case 8
			        SubFrame=[19];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[19];
			        RefFrame=[1];
			        SubJac=Jac_Global(1:3,:,ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_8(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_8(t,q,p,u,s,SubSubs);
			    case 9
			        SubFrame=[8];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-8];
			        RefFrame=[3];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_9(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_9(t,q,p,u,s,SubSubs);
			    case 10
			        SubFrame=[9];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-9];
			        RefFrame=[4];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_10(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_10(t,q,p,u,s,SubSubs);
			    case 11
			        SubFrame=[8];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-8];
			        RefFrame=[4];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_11(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_11(t,q,p,u,s,SubSubs);
			    case 12
			        SubFrame=[10];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-10];
			        RefFrame=[5];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_12(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_12(t,q,p,u,s,SubSubs);
			    case 13
			        SubFrame=[9];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-9];
			        RefFrame=[5];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_13(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_13(t,q,p,u,s,SubSubs);
			    case 14
			        SubFrame=[11];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-11];
			        RefFrame=[6];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_14(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_14(t,q,p,u,s,SubSubs);
			    case 15
			        SubFrame=[10];
			        SubSubs=0;
			        if SubFrame(1)~=0
			            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
			        end
			        ActFrame=[-10];
			        RefFrame=[6];
			        SubJac=Jac_Global(4:6,:,-ActFrame);
			        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_15(t,q,p,u,s,SubSubs);
			        SubGF=SubJac.'*SubEff;
			        SubInputJac=SubJac.'*InputEffJac_15(t,q,p,u,s,SubSubs);
			end
	        CollectGF(:,FCount)=SubGF;
	        CollectInputJac(:,:,FCount)=SubInputJac;
	    end
		function out1 = ForceJac_1(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_1
		%    OUT1 = FORCEJAC_1(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:22
		out1 = 0.0;
		end
		function out1 = ForceJac_2(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_2
		%    OUT1 = FORCEJAC_2(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:23
		out1 = 0.0;
		end
		function out1 = ForceJac_3(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_3
		%    OUT1 = FORCEJAC_3(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:25
		out1 = 0.0;
		end
		function out1 = ForceJac_4(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_4
		%    OUT1 = FORCEJAC_4(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:26
		out1 = 0.0;
		end
		function out1 = ForceJac_5(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_5
		%    OUT1 = FORCEJAC_5(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:28
		out1 = 0.0;
		end
		function out1 = ForceJac_6(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_6
		%    OUT1 = FORCEJAC_6(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:29
		out1 = 0.0;
		end
		function out1 = ForceJac_7(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_7
		%    OUT1 = FORCEJAC_7(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:31
		out1 = 0.0;
		end
		function out1 = ForceJac_8(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_8
		%    OUT1 = FORCEJAC_8(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:32
		out1 = 0.0;
		end
		function out1 = ForceJac_9(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_9
		%    OUT1 = FORCEJAC_9(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:34
		out1 = 0.0;
		end
		function out1 = ForceJac_10(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_10
		%    OUT1 = FORCEJAC_10(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:35
		out1 = 0.0;
		end
		function out1 = ForceJac_11(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_11
		%    OUT1 = FORCEJAC_11(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:37
		out1 = 0.0;
		end
		function out1 = ForceJac_12(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_12
		%    OUT1 = FORCEJAC_12(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:38
		out1 = 0.0;
		end
		function out1 = ForceJac_13(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_13
		%    OUT1 = FORCEJAC_13(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:40
		out1 = 0.0;
		end
		function out1 = ForceJac_14(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_14
		%    OUT1 = FORCEJAC_14(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:42
		out1 = 0.0;
		end
		function out1 = ForceJac_15(t,in2,in3,in4,NHSIGNAL,in6)
		%FORCEJAC_15
		%    OUT1 = FORCEJAC_15(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:43
		out1 = 0.0;
		end
		function out1 = InputEff_1(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_1
		%    OUT1 = INPUTEFF_1(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:22
		gAcc = in3(67,:);
		linkUkn3mpt1 = in4(5,:);
		out1 = [0.0;0.0;-gAcc.*linkUkn3mpt1];
		end
		function out1 = InputEff_2(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_2
		%    OUT1 = INPUTEFF_2(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:24
		gAcc = in3(67,:);
		linkUkn3mpt2 = in4(6,:);
		out1 = [0.0;0.0;-gAcc.*linkUkn3mpt2];
		end
		function out1 = InputEff_3(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_3
		%    OUT1 = INPUTEFF_3(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:25
		gAcc = in3(67,:);
		linkUkn3mpt3 = in4(7,:);
		out1 = [0.0;0.0;-gAcc.*linkUkn3mpt3];
		end
		function out1 = InputEff_4(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_4
		%    OUT1 = INPUTEFF_4(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:27
		gAcc = in3(67,:);
		linkUkn3mpt4 = in4(8,:);
		out1 = [0.0;0.0;-gAcc.*linkUkn3mpt4];
		end
		function out1 = InputEff_5(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_5
		%    OUT1 = INPUTEFF_5(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:28
		gAcc = in3(67,:);
		linkUkn4mpt1 = in4(9,:);
		out1 = [0.0;0.0;-gAcc.*linkUkn4mpt1];
		end
		function out1 = InputEff_6(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_6
		%    OUT1 = INPUTEFF_6(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:30
		gAcc = in3(67,:);
		linkUkn4mpt2 = in4(10,:);
		out1 = [0.0;0.0;-gAcc.*linkUkn4mpt2];
		end
		function out1 = InputEff_7(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_7
		%    OUT1 = INPUTEFF_7(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:31
		gAcc = in3(67,:);
		linkUkn4mpt3 = in4(11,:);
		out1 = [0.0;0.0;-gAcc.*linkUkn4mpt3];
		end
		function out1 = InputEff_8(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_8
		%    OUT1 = INPUTEFF_8(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:33
		gAcc = in3(67,:);
		linkUkn4mpt4 = in4(12,:);
		out1 = [0.0;0.0;-gAcc.*linkUkn4mpt4];
		end
		function out1 = InputEff_9(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_9
		%    OUT1 = INPUTEFF_9(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:34
		uJ1 = in4(1,:);
		out1 = [0.0;uJ1;0.0];
		end
		function out1 = InputEff_10(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_10
		%    OUT1 = INPUTEFF_10(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:36
		uJ2 = in4(2,:);
		out1 = [uJ2;0.0;0.0];
		end
		function out1 = InputEff_11(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_11
		%    OUT1 = INPUTEFF_11(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:37
		uJ2 = in4(2,:);
		out1 = [-uJ2;0.0;0.0];
		end
		function out1 = InputEff_12(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_12
		%    OUT1 = INPUTEFF_12(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:39
		uJ3 = in4(3,:);
		out1 = [0.0;uJ3;0.0];
		end
		function out1 = InputEff_13(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_13
		%    OUT1 = INPUTEFF_13(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:40
		uJ3 = in4(3,:);
		out1 = [0.0;-uJ3;0.0];
		end
		function out1 = InputEff_14(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_14
		%    OUT1 = INPUTEFF_14(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:42
		uJ4 = in4(4,:);
		out1 = [uJ4;0.0;0.0];
		end
		function out1 = InputEff_15(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFF_15
		%    OUT1 = INPUTEFF_15(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:44
		uJ4 = in4(4,:);
		out1 = [-uJ4;0.0;0.0];
		end
		function out1 = InputEffJac_1(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_1
		%    OUT1 = INPUTEFFJAC_1(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:23
		gAcc = in3(67,:);
		outSize=[3,12];
		elementRow=[3];
		elementCol=[5];
		elementList=[-gAcc];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_2(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_2
		%    OUT1 = INPUTEFFJAC_2(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:24
		gAcc = in3(67,:);
		outSize=[3,12];
		elementRow=[3];
		elementCol=[6];
		elementList=[-gAcc];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_3(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_3
		%    OUT1 = INPUTEFFJAC_3(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:26
		gAcc = in3(67,:);
		outSize=[3,12];
		elementRow=[3];
		elementCol=[7];
		elementList=[-gAcc];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_4(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_4
		%    OUT1 = INPUTEFFJAC_4(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:27
		gAcc = in3(67,:);
		outSize=[3,12];
		elementRow=[3];
		elementCol=[8];
		elementList=[-gAcc];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_5(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_5
		%    OUT1 = INPUTEFFJAC_5(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:29
		gAcc = in3(67,:);
		outSize=[3,12];
		elementRow=[3];
		elementCol=[9];
		elementList=[-gAcc];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_6(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_6
		%    OUT1 = INPUTEFFJAC_6(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:30
		gAcc = in3(67,:);
		outSize=[3,12];
		elementRow=[3];
		elementCol=[10];
		elementList=[-gAcc];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_7(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_7
		%    OUT1 = INPUTEFFJAC_7(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:32
		gAcc = in3(67,:);
		outSize=[3,12];
		elementRow=[3];
		elementCol=[11];
		elementList=[-gAcc];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_8(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_8
		%    OUT1 = INPUTEFFJAC_8(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:33
		gAcc = in3(67,:);
		outSize=[3,12];
		elementRow=[3];
		elementCol=[12];
		elementList=[-gAcc];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_9(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_9
		%    OUT1 = INPUTEFFJAC_9(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:35
		ang0x = in3(56,:);
		ang0y = in3(57,:);
		ang0z = in3(58,:);
		t2 = (cos(ang0z));
		t3 = (sqrt(3.0));
		t4 = (cos(ang0x));
		t5 = (sin(ang0x));
		t6 = (sin(ang0y));
		t7 = (sin(ang0z));
		t8 = (cos(ang0y));
		outSize=[3,12];
		elementRow=[1,2,3];
		elementCol=[1,1,1];
		elementList=[t3.*(t4.*t7-t2.*t5.*t6).*(-1.0./2.0)+(t2.*t8)./2.0,(t3.*(t2.*t4+t5.*t6.*t7))./2.0+(t7.*t8)./2.0,t6.*(-1.0./2.0)+(t3.*t5.*t8)./2.0];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_10(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_10
		%    OUT1 = INPUTEFFJAC_10(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:36
		ang0x = in3(56,:);
		ang0y = in3(57,:);
		ang0z = in3(58,:);
		qJ1__dt_0_ = in2(1,:);
		t2 = (sin(ang0x));
		t3 = (sin(ang0z));
		t4 = (cos(ang0x));
		t5 = (cos(ang0z));
		t6 = (sin(ang0y));
		t7 = (cos(qJ1__dt_0_));
		t8 = (sin(qJ1__dt_0_));
		t9 = (sqrt(3.0));
		t10 = (cos(ang0y));
		outSize=[3,12];
		elementRow=[1,2,3];
		elementCol=[2,2,2];
		elementList=[(t7.*(t3.*t4-t2.*t5.*t6))./2.0-t8.*(t2.*t3+t4.*t5.*t6)+(t5.*t7.*t9.*t10)./2.0,t7.*(t4.*t5+t2.*t3.*t6).*(-1.0./2.0)+t8.*(t2.*t5-t3.*t4.*t6)+(t3.*t7.*t9.*t10)./2.0,t2.*t7.*t10.*(-1.0./2.0)-t4.*t8.*t10-(t6.*t7.*t9)./2.0];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_11(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_11
		%    OUT1 = INPUTEFFJAC_11(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:38
		ang0x = in3(56,:);
		ang0y = in3(57,:);
		ang0z = in3(58,:);
		qJ1__dt_0_ = in2(1,:);
		t2 = (sin(ang0x));
		t3 = (sin(ang0z));
		t4 = (cos(ang0x));
		t5 = (cos(ang0z));
		t6 = (sin(ang0y));
		t7 = (cos(qJ1__dt_0_));
		t8 = (sin(qJ1__dt_0_));
		t9 = (sqrt(3.0));
		t10 = (cos(ang0y));
		outSize=[3,12];
		elementRow=[1,2,3];
		elementCol=[2,2,2];
		elementList=[t7.*(t3.*t4-t2.*t5.*t6).*(-1.0./2.0)+t8.*(t2.*t3+t4.*t5.*t6)-(t5.*t7.*t9.*t10)./2.0,(t7.*(t4.*t5+t2.*t3.*t6))./2.0-t8.*(t2.*t5-t3.*t4.*t6)-(t3.*t7.*t9.*t10)./2.0,(t2.*t7.*t10)./2.0+t4.*t8.*t10+(t6.*t7.*t9)./2.0];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_12(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_12
		%    OUT1 = INPUTEFFJAC_12(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:39
		ang0x = in3(56,:);
		ang0y = in3(57,:);
		ang0z = in3(58,:);
		qJ1__dt_0_ = in2(1,:);
		qJ2__dt_0_ = in2(2,:);
		t2 = (cos(ang0z));
		t3 = (cos(qJ2__dt_0_));
		t4 = (sqrt(3.0));
		t5 = (sin(qJ1__dt_0_));
		t6 = (sin(qJ2__dt_0_));
		t7 = (sin(ang0x));
		t8 = (sin(ang0z));
		t9 = (cos(ang0x));
		t10 = (sin(ang0y));
		t11 = ((t3.*t4)./2.0);
		t12 = (cos(ang0y));
		t13 = (t3./2.0);
		t14 = ((t4.*t5.*t6)./2.0);
		t15 = (t13+t14);
		t16 = (cos(qJ1__dt_0_));
		outSize=[3,12];
		elementRow=[1,2,3];
		elementCol=[3,3,3];
		elementList=[-(t11-(t5.*t6)./2.0).*(t8.*t9-t2.*t7.*t10)+t6.*t16.*(t7.*t8+t2.*t9.*t10)+t2.*t12.*t15,(t11-(t5.*t6)./2.0).*(t2.*t9+t7.*t8.*t10)-t6.*t16.*(t2.*t7-t8.*t9.*t10)+t8.*t12.*t15,-t10.*t15+t7.*t12.*(t11-(t5.*t6)./2.0)+t6.*t9.*t12.*t16];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_13(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_13
		%    OUT1 = INPUTEFFJAC_13(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:41
		ang0x = in3(56,:);
		ang0y = in3(57,:);
		ang0z = in3(58,:);
		qJ1__dt_0_ = in2(1,:);
		qJ2__dt_0_ = in2(2,:);
		t2 = (cos(ang0z));
		t3 = (cos(qJ2__dt_0_));
		t4 = (sqrt(3.0));
		t5 = (sin(qJ1__dt_0_));
		t6 = (sin(qJ2__dt_0_));
		t7 = (sin(ang0x));
		t8 = (sin(ang0z));
		t9 = (cos(ang0x));
		t10 = (sin(ang0y));
		t11 = ((t3.*t4)./2.0);
		t18 = ((t5.*t6)./2.0);
		t12 = (t11-t18);
		t13 = (cos(ang0y));
		t14 = (t3./2.0);
		t15 = ((t4.*t5.*t6)./2.0);
		t16 = (t14+t15);
		t17 = (cos(qJ1__dt_0_));
		outSize=[3,12];
		elementRow=[1,2,3];
		elementCol=[3,3,3];
		elementList=[t12.*(t8.*t9-t2.*t7.*t10)-t6.*t17.*(t7.*t8+t2.*t9.*t10)-t2.*t13.*t16,-t12.*(t2.*t9+t7.*t8.*t10)+t6.*t17.*(t2.*t7-t8.*t9.*t10)-t8.*t13.*t16,t10.*t16-t7.*t12.*t13-t6.*t9.*t13.*t17];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_14(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_14
		%    OUT1 = INPUTEFFJAC_14(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:43
		ang0x = in3(56,:);
		ang0y = in3(57,:);
		ang0z = in3(58,:);
		qJ1__dt_0_ = in2(1,:);
		qJ2__dt_0_ = in2(2,:);
		qJ3__dt_0_ = in2(3,:);
		t2 = (sin(qJ3__dt_0_));
		t3 = (sin(ang0x));
		t4 = (sin(ang0z));
		t5 = (cos(ang0x));
		t6 = (cos(ang0z));
		t7 = (sin(ang0y));
		t8 = (cos(qJ3__dt_0_));
		t9 = (sin(qJ1__dt_0_));
		t10 = (cos(qJ1__dt_0_));
		t11 = (cos(qJ2__dt_0_));
		t12 = (sin(qJ2__dt_0_));
		t13 = (sqrt(3.0));
		t14 = ((t2.*t12.*t13)./2.0);
		t15 = ((t2.*t9.*t11)./2.0);
		t25 = ((t8.*t10)./2.0);
		t16 = (t14+t15-t25);
		t17 = (t8.*t9);
		t18 = (t2.*t10.*t11);
		t19 = (t17+t18);
		t20 = (cos(ang0y));
		t21 = ((t2.*t12)./2.0);
		t22 = ((t8.*t10.*t13)./2.0);
		t24 = ((t2.*t9.*t11.*t13)./2.0);
		t23 = (t21+t22-t24);
		outSize=[3,12];
		elementRow=[1,2,3];
		elementCol=[4,4,4];
		elementList=[-t16.*(t4.*t5-t3.*t6.*t7)-t19.*(t3.*t4+t5.*t6.*t7)+t6.*t20.*t23,t16.*(t5.*t6+t3.*t4.*t7)+t19.*(t3.*t6-t4.*t5.*t7)+t4.*t20.*t23,-t7.*t23+t3.*t16.*t20-t5.*t19.*t20];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputEffJac_15(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTEFFJAC_15
		%    OUT1 = INPUTEFFJAC_15(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:44
		ang0x = in3(56,:);
		ang0y = in3(57,:);
		ang0z = in3(58,:);
		qJ1__dt_0_ = in2(1,:);
		qJ2__dt_0_ = in2(2,:);
		qJ3__dt_0_ = in2(3,:);
		t2 = (sin(qJ3__dt_0_));
		t3 = (sin(ang0x));
		t4 = (sin(ang0z));
		t5 = (cos(ang0x));
		t6 = (cos(ang0z));
		t7 = (sin(ang0y));
		t8 = (cos(qJ3__dt_0_));
		t9 = (sin(qJ1__dt_0_));
		t10 = (cos(qJ1__dt_0_));
		t11 = (cos(qJ2__dt_0_));
		t12 = (sin(qJ2__dt_0_));
		t13 = (sqrt(3.0));
		t14 = ((t2.*t12.*t13)./2.0);
		t15 = ((t2.*t9.*t11)./2.0);
		t25 = ((t8.*t10)./2.0);
		t16 = (t14+t15-t25);
		t17 = (t8.*t9);
		t18 = (t2.*t10.*t11);
		t19 = (t17+t18);
		t20 = (cos(ang0y));
		t21 = ((t2.*t12)./2.0);
		t22 = ((t8.*t10.*t13)./2.0);
		t24 = ((t2.*t9.*t11.*t13)./2.0);
		t23 = (t21+t22-t24);
		outSize=[3,12];
		elementRow=[1,2,3];
		elementCol=[4,4,4];
		elementList=[t16.*(t4.*t5-t3.*t6.*t7)+t19.*(t3.*t4+t5.*t6.*t7)-t6.*t20.*t23,-t16.*(t5.*t6+t3.*t4.*t7)-t19.*(t3.*t6-t4.*t5.*t7)-t4.*t20.*t23,t7.*t23-t3.*t16.*t20+t5.*t19.*t20];
		out1=zeros(outSize(1),outSize(2));
		for eleNum=1:length(elementList)
		    out1(elementRow(eleNum),elementCol(eleNum))=elementList(eleNum);
		end
		end
		function out1 = InputFVJac_1(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_1
		%    OUT1 = INPUTFVJAC_1(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:23
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_2(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_2
		%    OUT1 = INPUTFVJAC_2(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:24
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_3(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_3
		%    OUT1 = INPUTFVJAC_3(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:26
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_4(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_4
		%    OUT1 = INPUTFVJAC_4(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:28
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_5(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_5
		%    OUT1 = INPUTFVJAC_5(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:29
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_6(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_6
		%    OUT1 = INPUTFVJAC_6(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:30
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_7(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_7
		%    OUT1 = INPUTFVJAC_7(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:32
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_8(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_8
		%    OUT1 = INPUTFVJAC_8(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:33
		out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
		end
		function out1 = InputFVJac_9(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_9
		%    OUT1 = INPUTFVJAC_9(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:35
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
		function out1 = InputFVJac_10(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_10
		%    OUT1 = INPUTFVJAC_10(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:36
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
		function out1 = InputFVJac_11(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_11
		%    OUT1 = INPUTFVJAC_11(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:38
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
		function out1 = InputFVJac_12(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_12
		%    OUT1 = INPUTFVJAC_12(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:40
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
		function out1 = InputFVJac_13(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_13
		%    OUT1 = INPUTFVJAC_13(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:41
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
		function out1 = InputFVJac_14(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_14
		%    OUT1 = INPUTFVJAC_14(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:43
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
		function out1 = InputFVJac_15(t,in2,in3,in4,NHSIGNAL,in6)
		%INPUTFVJAC_15
		%    OUT1 = INPUTFVJAC_15(T,IN2,IN3,IN4,NHSIGNAL,IN6)
		%    This function was generated by the Symbolic Math Toolbox version 8.2.
		%    07-Jan-2021 01:58:45
		out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
		end
	end

	function [ConsJac,ConsCor,ConsGF]=Constraint_Arm4D(t,q,p,u,s,TransDis_Global,Vel_Global,Cor_Global,Jac_Global,Quat_Global)
	%% Constraint dynamic property calculator (Toolbox Internal Use)
	%% 
	%% Robotics & Mechatronics Lab, Virginia Tech. (No Copyright Claimed)
	%% Author: Jiamin Wang; Revised: 28-Jan-2019
	    ConsNum=0;
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
	        
	        ConsJac(ConsCount,:)=SubConsJac;
	        ConsCor(ConsCount,:)=SubConsCor;
	        ConsGF(ConsCount,:)=SubConsGF;
	    end
	end

end

