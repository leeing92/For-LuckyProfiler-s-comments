function [Result,SingleResult,FinalTro10Id,PreProjectionNum]=How2FWHM2Auto(RenderedImg,px,SaveLabel,Method)
%% This code is to calculate the full image and Point where has the reliable resolution belong several FWHM method
%Input: 
%   RenderedIm--Img;
%   px--pixle size
%Output:
%   SingelResult:[SkelID,Locx,Locy,angle,rr,SingleProfile]
%   Result:{SkelID,Loc,ProjectionProfile,CalculatedResolution,SolvedWidth,Fitresult,Fitgood}
%% Prepare for whole image 
%     str=['E:\Mine\LAB\Resolution\FWHM method for LM\Auto-Projection\Other structure\',SaveLabel];
%     str=['E:\Mine\LAB\Resolution\FWHM method for LM\WhyNeed\Supplement For Paper\Density\',SaveLabel];
% str=['E:\Mine\LAB\Resolution\FWHM method for LM\Auto-Projection\Auto-CadidateSelect-20220405\',SaveLabel];
    if size(RenderedImg,3)==3
        RenderedImg=rgb2gray(RenderedImg);
    end
    SizeIm1=size(RenderedImg,1);
    SizeIm2=size(RenderedImg,2);
    Offset=[RenderedImg(1,:),RenderedImg(SizeIm1,:),RenderedImg(:,1)',RenderedImg(:,SizeIm2)'];
    II=RenderedImg-mean(Offset,2);
    n=ceil(20/px);
    se=strel('disk',n);
    Im=DenoiseFunc(II,px); 
    RenderedImg2BS=imdilate(Im,se); 
    Bounderpre=imbinarize(RenderedImg2BS,'adaptive','Sensitivity',0.4); 
    h=fspecial('laplacian',0);
    Im2select=imfilter(Bounderpre,h,'same');
    se=strel('disk',1);
    Im2select=imdilate(Im2select,se); 
    thres=30;
    skeletontemp1=Normalized(Im);
%     skeletontemp=regionThres(Im);
    skeletontemp=imbinarize(skeletontemp1,'adaptive','Sensitivity',0.4);
    [skr,rad] = skeleton( skeletontemp);
    SkelImwithout1 = bwmorph(skr > thres,'skel',inf);
    SkelNum=sum(SkelImwithout1(:));
    [Locx,Locy]=find(SkelImwithout1);
    

%% Calcualte each location to get r and lineprofile    
    [LineProfile,ang,rr]=EachProfile(Locx,Locy,SkelNum,SkelImwithout1,Im2select,II);
    rrmax=max(rr(:))+1;
    %start-add on 20220405
    rrnozero=rr(rr~=0);
    [a1,a2]=histcounts(rrnozero);
    a3=(a2(2:end)+a2(1:end-1))/2;
    fitresult=fit(a3',a1'/sum(a1),'Gauss1');
    T1=fitresult.b1+3*fitresult.c1;
    T2=fitresult.b1-3*fitresult.c1;
    Index1=rr>T1 | rr<T2;
    rr(Index1)=0;
    %end-add on 20220405
    
     SingleResult=zeros(SkelNum,2*rrmax+1+5);
     for j=1:SkelNum         
        SingleResult(j,1)=j;
        index1=Locx(j);
        index2=Locy(j);
        SingleResult(j,2:3)=[index1,index2];
        SingleResult(j,4)=ang(j);
        SingleResult(j,5)=rr(j);
        if rr(j)>0
        rradd=rrmax-rr(j);
        SingleResult(j,6:size(SingleResult,2))=[zeros(1,rradd),LineProfile{j}(:)',zeros(1,rradd)];
 %test---------------
%     figure(1);
%     plot(SingleResult(j,4:size(SingleResult,2)));hold on
%test end----------------
        end
     end
%     hold off    
%% According r selecte ROI      
     %Pre-thickness setting (Projection Num)
     if SkelNum>100
         randnum=100;
     else
         randnum=SkelNum;
     end
     [PreProjectionNum,intenum]=PreProNum(randnum,SkelImwithout1,SkelNum,SingleResult,px,SizeIm1,SizeIm2);
    
    % pre-r for each Region
    parfor i=1:SkelNum  
        [RCadi(i,:),RelatedLocIndex(i,:)]=PreR(SkelImwithout1,SingleResult,rr,i,PreProjectionNum,SizeIm1,SizeIm2);       
    end
    
    % According r choosing and calculate resolution
    %start-add on 20220405
    for i=1:size(RCadi,1)
        RCadiTemp=RCadi(i,:);
        if sum(RCadiTemp~=0)
        RCadiMean(i)=mean(RCadiTemp(RCadiTemp~=0));
        else
            RCadiMean(i)=0;
        end
    end
    [R,ind]=sort(RCadiMean);
 %end-add on 20220405
 
%     RCadiMean=mean(RCadi,2);
%     [R,ind]=sort(RCadiMean);
    Index=find(R~=0);
    N=ceil(0.04*SkelNum); %changed by 20220407
    SkelId=ind(Index(1):Index(N));
    for i=1:N
        SI=SkelId(i);
        LocIndex=RelatedLocIndex(SI,:);
        for j=1:size(LocIndex,2)
            LP(j,:)=SingleResult(LocIndex(j),6:size(SingleResult,2));
         end
        L=Normalized(sum(LP));
        switch Method
            case 'Size'
            [WithWidthFit_s,WithWidthFit_r,h]=WidthSolvedFWHM1D(L',px);
            Result{i,1}=SI;
            Result{i,2}=[LocIndex',SingleResult(LocIndex,2:3)];
            Result{i,3}=L;
            Result{i,4}=WithWidthFit_s;
            Result{i,5}=WithWidthFit_r;
            case 'Gauss'
                Lsize=px*(1:length(L));
            [fitresult,gof]=fit(Lsize',L','Gauss1','Lower',[1 0 0],'Upper',[1 Inf Inf]);
            Result{i,1}=SI;
            Result{i,2}=[LocIndex',SingleResult(LocIndex,2:3)];
            Result{i,3}=L;
            Result{i,4}=fitresult.c1*1.665;
            Result{i,5}=gof.rsquare;
%             h=figure;
%             plot(fitresult,Lsize',L');
        end
     %  Show Line Profile
    %-----------
%     maxL=max(LP(:));
%     Lsize=1:length(L);
%     Lsize=Lsize*px;
%     figure(h);
%     for jj=1:10
%          plot(Lsize,LP(jj,:)/maxL,':');
%     end
%     plot(Lsize,L,'Color','k','MarkerSize',3,'LineWidth',2);
%     ylabel('Normalized Intensity','FontSize',16,'FontWeight','bold');
%     xlabel('Distance (nm)','FontSize',16,'FontWeight','bold');
%     set(gca,'FontSize',16,'LineWidth',2); 
%     box off;
%     hold off
%     savefig([str,'-',num2str(i),'.fig']);
    %----------
    %  Show Choosing Loc
%     Im2Show=ShowROI(Result,SizeIm1,SizeIm2);
%     figure;imshow(Im2Show);
%     savefig([str,'-Loc',num2str(i),'.fig']);
    close all
    end
    [FinalTro10Id]=ChoosOne(Result);
    disp(['Final ID Through ',num2str(N),' Candidate is ',num2str(FinalTro10Id)]);
    Im2Show=ShowROI(Result(FinalTro10Id,:),SizeIm1,SizeIm2);
    figure;imshow(Im2Show);
    Lfinal=Result{FinalTro10Id,3};
    xlength=1:length(Lfinal);
    xlength=xlength*px;
    [WithWidthFit_s,WithWidthFit_r,h]=WidthSolvedFWHM1D(Lfinal',px);
    fitresult=fit(xlength',Lfinal','Gauss1','Lower',[1 0 0],'Upper',[1 Inf Inf]);
    figure(h);
    plot(fitresult,xlength',Lfinal);
    ylabel('Normalized Intensity','FontSize',16,'FontWeight','bold');
    xlabel('Distance (nm)','FontSize',16,'FontWeight','bold');
    set(gca,'FontSize',16,'LineWidth',2); 
    box off;
    legend('Deconvolution','Raw','Gaussian');
%     save([str,'.mat']);
end

function [LineProfile,ang,rr]=EachProfile(Locx,Locy,SkelNum,SkelImwithout1,Im2select,II)
    rr=zeros(1,SkelNum);
    LineProfile=cell(1,SkelNum);
    ang=zeros(1,SkelNum);
%     check=[279;284;294;299;308;313;330;345;368;385;391;403;409;424;439;621;665;770;789;803;832;920;1179;1195;1196;1212;1277;1291;1447;1451;1466;1477;1481;1486;1492;1496;1509;1512;1523;1524;1527;1544;1545;1547;1561;1562;1565;1577;1578;1581;1593;1594;1610;1626;1645;1665;1683;1700;1708;1716;1718;1734;1750;1757;1766;1774;1783;1790;1800;1809;1819;1826;1830;1844;1849;1884;1891;1902;1909;1916;1925;1931;1932;1940;1948;1949;1951;1957;1963;1964;1978;1980;1997;1999;2014;2017;2023;2029;2032;2045;2048;2061;2064;2079;2085;2095;2101;2117;2135;2152;2162;2168;2177;2185;2194;2203;2213;2221;2230;2247;2265;2282;2297;2298;2316;2317;2324;2333;2341;2351;2360;2368;2376;2385;2386;2393;2404;2420;2437;2455;2462;2479;2499;2515;2524;2531;2541;2548;2565;2611;2628;2635;2646;2653;2664;2681;2697;2704;2713;2728;2731;2748;2756;2767;2784;2797;2801;2803;2809;2814;2817;2825;2826;2831;2834;2843;2847;2850;2858;2862;2866;2874;2878;2882;2891;2897;2905;2912;2920;2921;2928;2938;2955;2973;2991;2998;3007;3024;3038;3042;3060;3066;3077;3085;3095;3096;3102;3111;3118;3126;3134;3142;3149;3158;3165;3180;3194;3202;3211;3218;3226;3233;3234;3241;3248;3256;3271;3310;3312;3313;3325;3328;3335;3339;3342;3348;3355;3365;3374;3380;3394;3416;3431;3439;3446;3463;3523;3539;3554;3555;3558;3562;3563;3570;3571;3574;3586;3587;3590;3601;3605;3615;3620;3659;3672;3675;3687;3688;3705;3724;3756;3772;3819;3835;3851;3867;3883;3884;3971;3990;4006;4007;4023;4024;4039;4040;4055;4056;4071;4072;4086;4087;4104;4105;4120;4137;4207;4224;4241;4260;4261;4277;4278;4293;4302;4312;4330;4372;4383;4400;4420;4438;4439;4444;4458;4481;4498;4499;4517;4537;4557;4618;4635;4636;4656;4657;4675;4714;4735;4753;4754;4774;4775;4785;4795;4800;4818;4837;4857;4871;4876;4890;4895;4901;4909;4913;4920;4932;4938;4955;4959;4967;4973;4992;5012;5032;5037;5038;5046;5050;5055;5056;5063;5064;5073;5074;5081;5082;5091;5100;5109;5117;5130;5141;5148;5161;5162;5169;5181;5182;5191;5199;5200;5202;5217;5220;5236;5240;5256;5259;5275;5277;5284;5293;5295;5314;5335;5354;5372;5391;5393;5396;5409;5428;5432;5450;5469;5487;5502;5504;5521;5538;5557;5559;5575;5593;5611;5649;5662;5667;5683;5701;5713;5716;5717;5731;5734;5735;5750;5751;5769;5786;5838;5850;5854;5866;5870;5884;5901;5920;5936;5937;5941;5959;5978;5991;6008;6026;6048;6065;6078;6083;6087;6097;6102;6120;6124;6138;6142;6155;6171;6191;6208;6221;6225;6242;6332;6345;6348;6366;6382;6400;6412;6413;6416;6431;6434;6447;6450;6464;6467;6481;6499;6515;6533;6549;6564;6579;6596;6611;6626;6642;6659;6678;6736;6740;6753;6757;6771;6775;6793;6814;6831;6846;6850;6865;6869;6884;6888;6900;6904;6917;6921;6937;6941;6959;6974;6993;7010;7022;7025;7037;7040;7102;7118;7133;7150;7165;7181;7198;7267;7426;7442;7459;7476;7494;7511;7529;7546;7563;7570;7583;7588;7600;7623;7635;7642;7661;7679;7801;7811;7813;7831;7863;7882;8061;8161;8179;8304;8484;8498;8503;8519;8520;8540;8724;8743;8762;8763;8782;9089;9259;9273;9316;9335;9347;9352;9369;9419;9436;9452;9468;9682;9702;9719;9737;9738;9751;9821;9831;9948;9966;9983;10001;10019;10037;10042;10054;10059;10060;10071;10076;10087;10104;10122;10140;10228;10229;10246;10247;10263;10264;10281;10298;10316;10333;10351;10369;10386;10387;10405;10424;10441;10483;10511;10527;10532;10545;10670;10685;10702;10718;10738;10756;10794;10812;10830;10847;10865;10918;10935;10951;10968;10984;11001;11045;11061;11078;11083;11102;11118;11119;11136;11147;11163;11184;11191;11201;11208;11217;11228;11233;11249;11272;11303;11336;11344;11353;11367;11369;11376;11379;11380;11396;11397;11398;11410;11411;11424;11437;11451;11465;11480;11494;11498;11509;11512;11514;11524;11529;11530;11540;11544;11545;11555;11560;11561;11572;11577;11581;11587;11596;11602;11611;11617;11622;11633;11638;11648;11663;11720;11730;11734;11745;11749;11759;11763;11773;11777;11788;11792;11803;11807;11816;11818;11822;11831;11836;11845;11851;11868;11883;11901;11917;11931;11934;11947;11950;11967;11977;11982;11998;12009;12012;12023;12026;12041;12057;12074;12087;12088;12091;12104;12116;12122;12137;12153;12168;12186;12188;12189;12193;12205;12220;12237;12238;12242;12253;12257;12272;12292;12293;12310;12311;12313;12335;12353;12368;12372;12386;12390;12404;12408;12422;12426;12439;12441;12460;12475;12477;12519;12554;12576;12596;12597;12601;12614;12615;12635;12654;12657;12670;12691;12708;12721;12725;12726;12738;12742;12743;12759;12776;12781;12794;12795;12811;12812;12815;12828;12829;12845;12846;12863;12883;12899;12915;12916;12917;12930;12933;12934;12951;13015;13023;13027;13040;13041;13048;13057;13060;13079;13098;13134;13150;13179;13197;13215;13272;13295;13296;13310;13339;13354;13368;13374;13383;13389;13397;13414;13430;13443;13507;13522;13537;13552;13592;13607;13623;13628;13639;13644;13654;13660;13670;13675;13685;13691;13702;13707;13711;13717;13726;13732;13741;13747;13756;13757;13763;13771;13772;13778;13783;13794;13799;13808;13813;13823;13828;13876;13898;13904;13916;13919;13922;13937;13942;13949;13953;13954;13957;13964;13968;13969;13973;13980;13984;13985;13989;13996;14000;14001;14004;14015;14016;14020;14031;14032;14035;14047;14048;14052;14063;14064;14067;14078;14083;14095;14099;14109;14113;14123;14127;14138;14153;14167;14176;14184;14194;14203;14222;14237;14242;14247;14257;14262;14270;14272;14279;14285;14286;14288;14299;14301;14309;14316;14324;14330;14338;14343;14351;14357;14365;14370;14379;14385;14398;14400;14404;14409;14417;14425;14426;14432;14440;14447;14454;14467;14485;14536;14603;14674;14712;14879;14953;14958;14970;14976;15022;15028;15154;15166;15168;15376;15560;15575;15590;15624;15673;15674;15690;15694;15695;15707;15723;15740;15882;15900;15907;15911;15989;15993;16007;16012;16027;16031;16045;16050;16062;16066;16078;16082;16095;16099;16103;16113;16118;16122;16130;16134;16138;16147;16151;16155;16163;16169;16173;16183;16187;16191;16206;16256;16260;16264;16273;16276;16278;16293;16310;16324;16326;16328;16333;16340;16342;16361;16363;16379;16407;16424;16441;16482;16499;16506;16512;16513;16515;16522;16528;16532;16539;16544;16547;16554;16561;16572;16578;16589;16595;16597;16605;16610;16612;16621;16627;16629;16637;16642;16644;16652;16658;16668;16673;16684;16690;16702;16709;16719;16725;16735;16742;16752;16758;16769;16784;16825;16906;16918;16919;16921;16923;16933;16934;16936;16938;16950;16952;16979;17020;17041;17044;17060;17096;17122;17134;17136;17147;17149;17202;17216;17230;17233;17316;17403;17419;17432;17447;17511;17526;17543;17545;17558;17577;17578;17592;17607;17655;17755;17771;17788;17804;17821;17837;17856;17874;17890;17907;17923;17972;17990;18007;18011;18091;18315;18330;18521;18543;18565;18580;18595;18609;18625;18642;18658;18675;18688;18703;18720;18735;18750;18751;18752;18767;18785;18787;18802;18818;18834;18850;18866;18884;18892;18900;18919;18935;18952;18954;18967;18969;18982;18984;18989;18998;19001;19007;19014;19021;19031;19037;19038;19049;19055;19056;19075;19103;19105;19125;19141;19152;19158;19167;19175;19181;19190;19197;19220;19237;19253;19269;19287;19306;19321;19371;19385;19427;19455;19468;19528];
    for i= 1:SkelNum % check(1):check(length(check))
%        below judges angle through 8 ajdent pixels
       [angle,G]=adjangle1(Locx,Locy,i,SkelImwithout1);
       ang(i)=angle;
       BoundCount=0;
       LineProfile_r=[];
       LineProfile_l=[];
       XX1=0;YY1=0;XX2=0;YY2=0;rep=0;
       for r=1:1000
            X1=Locx(i)-r*sin(angle); %suit with new angle
           Y1=Locy(i)+r*cos(angle);
           X1=round(X1);
           Y1=round(Y1);
           X2=Locx(i)-r*sin(pi+angle); %suit with new angle
           Y2=Locy(i)+r*cos(pi+angle);
           X2=round(X2);
           Y2=round(Y2);
           index1=Locx(i);
           index2=Locy(i);
           if Y1>size(Im2select,2)|| Y1<1 || X1>size(Im2select,1)||X1<1||Y2>size(Im2select,2)|| Y2<1 || X2>size(Im2select,1)||X2<1
%                rr(i)=r;
%                 [Locx(i),Locy(i);X1,Y1;X2,Y2]
               break
           end
           if (X1==XX1 && Y1==YY1) || (X2==XX2 && Y2==YY2)
               rep=rep+1;
           else
                XX1=X1;YY1=Y1;XX2=X2;YY2=Y2;
                BoundCount=Im2select(X1,Y1)+Im2select(X2,Y2)+BoundCount;
                LineProfile_r=[LineProfile_r,II(X1,Y1)];
                LineProfile_l=[II(X2,Y2),LineProfile_l];
               if BoundCount>=2
                rr(i)=r-rep;
                    break
               end
           end
       end
       LineProfile{i}=[LineProfile_l,II(index1,index2),LineProfile_r];
%        disp(['SkelNum:',num2str(i)]);
%         [Locx(i),Locy(i);X1,Y1;X2,Y2]
%        figure(3);plot(LineProfile{i});
    end
end

function [PreProjectionNum,intenum]=PreProNum(randnum,SkelImwithout1,SkelNum,SingleResult,px,SizeIm1,SizeIm2)
     PreIndex1=round(SkelNum*rand(1,randnum));
     aTry=SingleResult(:,6:size(SingleResult,2));
     Cannot=find(sum(aTry,2)==0);
     Cannum=0;
     for i=1:randnum
         if sum(PreIndex1(i)==Cannot)
             Cannum=Cannum+1;
         else
             PreIndex(i-Cannum)=PreIndex1(i);
         end
     end
     randnum=length(PreIndex);
     intenum=zeros(1,randnum);
    parfor i=1:randnum
        Index=PreIndex(i);
        XLoc=SingleResult(Index,2);
        YLoc=SingleResult(Index,3);
        Loc=[XLoc,YLoc];
        prefitg=0;
        LP=[];
        L=SingleResult(Index,6:size(SingleResult,2));
        L=Normalized(L);
        Lsize=1:length(L);
        Lsize=Lsize*px;
%         figure(1);plot(Lsize,L);
        [fitresult,fitgood]=fit(Lsize',L','Gauss1');        
        prefitg=fitresult.c1*1.665;
      if XLoc-1>=1 && YLoc-1>=1 && XLoc+1<=SizeIm1 && YLoc+1<=SizeIm2 
        for k=1:1000    
            Loc=aggregate(Loc,SkelImwithout1);          
            for j=1:size(Loc,1)
                Ind=find(SingleResult(:,2)==Loc(j,1) & SingleResult(:,3)==Loc(j,2));
                LP(j,:)=SingleResult(Ind,6:size(SingleResult,2));
            end
            L=Normalized(sum(LP,1));
            [fitresult,fitgood]=fit(Lsize',L','Gauss1');
            fg=fitgood.rsquare;
            sig=fitresult.c1*1.665;
            if fg>0.8 && abs(sig-prefitg)<0.01
                intenum(i)=size(Loc,1);
                break
            else
                prefitg=sig;
            end
        end
      end
    end
    PreProjectionNum=ceil(sum(intenum)/randnum); % changed
%     [counts,centers] = hist(intenum);
%     [~,ind]=sort(counts,'descend');
%     TempId=find(ind==ind(2));
%     TempId=[ind(1),TempId];
%     PreProjectionNum=ceil(mean(centers(TempId)));
end

function  [RCadi,Index]=PreR(SkelImwithout1,SingleResult,rr,i,PreProjectionNum,SizeIm1,SizeIm2)
        XLoc=SingleResult(i,2);
        YLoc=SingleResult(i,3);
        Loc=[XLoc,YLoc];
      if XLoc-1>=1 && YLoc-1>=1 && XLoc+1<=SizeIm1 && YLoc+1<=SizeIm2 
        for k=1:1000    
            Loc=aggregate(Loc,SkelImwithout1);
            if size(Loc,1)>=PreProjectionNum
%                 jmax=size(Loc,1);
%                 for j=1:jmax
                 for j=1:PreProjectionNum
                    Index(j)=find(SingleResult(:,2)==Loc(j,1)&SingleResult(:,3)==Loc(j,2));
                    RCadi(j)=rr(Index(j));
                end
                break
            else
                Index=zeros(1,PreProjectionNum);
                RCadi=zeros(1,PreProjectionNum);
            end
        end
       else
          Index=zeros(1,PreProjectionNum);
          RCadi=zeros(1,PreProjectionNum);
      end
end

function [FinalTro10Id]=ChoosOne(Result)
    A=Result(:,4);
    A=cell2mat(A);
%     PreSort=abs(mean(A)-A); % changed
%     [~,Id]=sort(PreSort);
    [~,Id]=sort(A);
    FinalTro10Id=Id(1);
end
