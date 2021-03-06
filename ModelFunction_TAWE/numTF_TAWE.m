function [frameTF]=numTF_TAWE(t,q,p,u,s)
    TFcalc=linkTF(t,q,p,u,s);
    FrameNum=26;
    TF=reshape(TFcalc,[4,4,FrameNum]);
    frameTF=zeros(4,4,FrameNum);
    frameTF(:,:,1)=TF(:,:,1);
    frameCheckList=zeros(1,FrameNum);
    frameCheckList(1)=1;
    for frame=2:FrameNum
        framePath=0;
        %thisTF=eye(4,4);
        
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

        [frameCheckList,TF,frameTF]=frameTFRecur(framePath,frameCheckList,TF,frameTF);
        %{
        for jj=numel(framePath):-1:1
            curNum=framePath(jj);
            curTF=TF(:,:,curNum);
            thisTF=curTF*thisTF;
        end
        
        frameTF(:,:,frame)=thisTF;
        %}
    end
    function [frameCheckList_,TF_,frameTF_]=frameTFRecur(framePath_,frameCheckList_,TF_,frameTF_)
        curNum_=framePath_(end);
        if(~frameCheckList_(curNum_))
            preNum_=framePath_(end-1);
            if(frameCheckList_(preNum_))
                frameTF_(:,:,curNum_)=frameTF_(:,:,preNum_)*TF_(:,:,curNum_);
                frameCheckList_(curNum_)=1;
            else
                [frameCheckList_,TF_,frameTF_]=frameTFRecur(framePath_(1:end-1),frameCheckList_,TF_,frameTF_);
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

end

