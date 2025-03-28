%%%Input data:(1) image of speckle, 
% (2) binary image of pol2 location (a single pixel selected for each pol2 using ImageJ "Find Maxima" and generate binary image). 

numbofimg=10; %modify if you have consecutive number labeled for your image.
for N=1:numbofimg %
    NN=num2str(N);
    %Modify directory, file name and format
    dir='/directory for images/'; 
    fname=['your file name' NN  ];format='.jpg';
    fname2=['your file name' NN ]; format2='.tif';

    %%Thrap3 segmentation 'Find good value for "thr"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Th1=imread([dir fname format]);
    %Th1=imbinarize(Th1); % for global thresholding
    T = adaptthresh(I,thr);
    BW = imbinarize(I,T);

    %%Ser2P bw image%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %S2P=imbinarize(imread([dir fname2 format2]));
    S2P=imread([dir fname2 format2]);
    S2P=ones(size(S2P))-imbinarize(S2P); %Ser2P image has 1 white in background.
    %figure, imshow(I2)
    %%%%Thrap3=3 Ser2P=2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s0=Th1*3+S2P*2;
    sz=size(s0)+10; %line, column (227,264)
    s=zeros(sz);
    s(5:sz(1,1)-6, 5:sz(1,2)-6)=s0;
    s0large=s;
    %%%%%%1> Count Ser2P on Thrap3, s+t=5
    f=find(s==5); %%% need to add in the end 
    S2P_ovlp=length(f); %Ser2P on Thrap3
    %turn the Ser2P+Thrap3 into Thrap3 value.
    s(f)=3;

    %%%%%%2> Ser2P distance from Thrap3
    %To avoid errors, add frame sized r0 =size of mask on image.
    r0 =5;%7;%5;%radius of mask
    f1=find(s==2); %location of single Ser2P foci
 
    %%%Lets make a mask %%size--> #line x #column
    %coordinate (x)= Column
    xc =ceil(f1/sz(1,1));
    % floor ()+1
    %coordinate (y)= Line
    yc =(f1/sz(1,1)-floor(f1/sz(1,1)))*sz(1,1);% f1(:,2); %edit from here...
    x_ct3=zeros(length(f1),1);
    y_ct3=zeros(length(f1),1);

    for fi=1:length(f1)
        fi=fi
        [xx yy] = meshgrid(1:sz(2), 1:sz(1));%for tif %meshgrid(1:sI_r(1), 1:sI_r(2));
        % Next create the circle in the image.
        mask0= (xx - xc(fi)).^2 + (yy - yc(fi)).^2 <= r0^2;%Imagej x,y coord=> matlab y,x  
        mask0i=imcomplement(mask0); %reverse B&W 
        surr_area=mask0i*99999;%double
        Im=s-surr_area;
        Im(find(Im<0))=0;%s<0))=0;
        
        %Location of S2P
        %xc(fi), yc(fi) %find (Im==2)
        %Location of Th3
        loc_t3=find(Im==3);
        
        if isempty(loc_t3)==1 
            r1=20;
            mask0= (xx - xc(fi)).^2 + (yy - yc(fi)).^2 <= r1^2;%Imagej x,y coord=> matlab y,x  
            mask0i=imcomplement(mask0); %reverse B&W 
            surr_area=mask0i*99999;%double
            Im=s-surr_area;
            Im(find(s<0))=0;
            %figure, imshow(Im)
            
        end    
        
         loc_t3=find(Im==3);
         xc_t3 =ceil(loc_t3/sz(1,1));
         yc_t3 =(loc_t3/sz(1,1)-floor(loc_t3/sz(1,1)))*sz(1,1);
      %find the shortest path
        dist=zeros(length(loc_t3),1);
        for t=1:length(loc_t3)
            dist(t,1)=((xc(fi)-xc_t3(t))^2+(yc(fi)-yc_t3(t))^2)^0.5;
        end
        mdist=min(dist);
        dpx(fi,1)=min(dist);
        
        loc_dpx=find(dist==dpx(fi));
            if length (loc_dpx)>1
                loc_dpx=1;
            end    
            
        x_ct3(fi,1)=xc_t3(loc_dpx);
        y_ct3(fi,1)=yc_t3(loc_dpx);
%         figure, imshow(Im)
%         RGB = insertMarker(s,[xc(fi) ,yc(fi) ],'s','color','red','size', 1);
%         RGB1 = insertMarker(RGB,[x_ct3 ,y_ct3 ],'s','color','green','size', 1);
%         figure, imshow(RGB1)  
   
    end
% % %final check    
%         RGBf = insertMarker(s,[xc, yc],'s','color','red','size', 1);
%         RGB1f = insertMarker(RGBf,[x_ct3, y_ct3],'s','color','green','size', 1);
%         figure, imshow(RGB1f)
%         hold on, title ('S2P with d=0 not shown')

%final data
    dist_btwS2P_T3=dpx*0.079; %pixel size = 0.079um
    s2p_onT3=S2P_ovlp;
    d_ovlp=zeros(s2p_onT3,1);
    

%%% Lets show pol2 in different color depending on distance..
       %%%[1] x, y coordinate for pol2 on speckle
    s0s=size(s0large);
    xc0 =ceil(f/s0s(1,1));
    yc0 =(f/s0s(1,1)-floor(f/s0s(1,1)))*s0s(1,1);    
    total_x=[xc;xc0];
    total_y=[yc;yc0];
    total_distances=[dist_btwS2P_T3;d_ovlp];
    %figure, histogram(total_distances)
    
    numb_of_x=length(total_x)
    numb_of_y=length(total_y) 
    numb_of_orgS2P=length(find(S2P==1))
    
    numb_of_dist=length(total_distances)
    
   Ass= find(total_distances<0.15);
   Notsure= find(total_distances==0.15 | total_distances>0.15 & total_distances<0.45);
   NA= find(total_distances==0.45 | total_distances>0.45);
    
    A_x=total_x(Ass);
    Notsure_x=total_x(Notsure);
    NA_x=total_x(NA);
    
    A_y=total_y(Ass);
    Notsure_y=total_y(Notsure);
    NA_y=total_y(NA);

    % figure, imshow(Im)

    %Save final data (collected distance measurements)   
    fname_text=[dir fname2 '_dist.txt'];
    fid = fopen(fname_text,'w');
    fprintf(fid, '%12.8f\n', total_distances);
    fclose(fid);
    
    clear
end



