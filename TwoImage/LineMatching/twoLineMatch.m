function [linematch1, linematch2] = twoLineMatch(imgpath1, imgpath2, pts1, pts2, parameters) 

img1 = imread(imgpath1);  img2 = imread(imgpath2);

lines_1 = lsd(imgpath1);  lines_1 = refineLine(lines_1, img1);
lines_2 = lsd(imgpath2);  lines_2 = refineLine(lines_2, img2);% detect line segments
len_lines1 = sqrt( (lines_1(2,:)-lines_1(1,:)).^2 + (lines_1(4,:)-lines_1(3,:)).^2 );
len_lines2 = sqrt( (lines_2(2,:)-lines_2(1,:)).^2 + (lines_2(4,:)-lines_2(3,:)).^2 ); % length of line segments


% choose the first 200 line segments
%length_t = min([500, length(len_lines1), length(len_lines2)]); 
len_threshold = parameters.line_threshold;

lines1 = lines_1(1:4,len_lines1>=len_threshold);  % target
lines2 = lines_2(1:4,len_lines2>=len_threshold);  % reference

lines1 = [lines1(1,:)',lines1(3,:)', lines1(2,:)', lines1(4,:)'];
lines2 = [lines2(1,:)',lines2(3,:)', lines2(2,:)', lines2(4,:)'];

disp(' Reading files and preparing...');
[lines1, pointlist1]=paras(img1, lines1);
[lines2, pointlist2]=paras(img2, lines2);

    len1=length(lines1);
    len2=length(lines2);
    sublinds1=1:length(lines1);
    sublinds2=1:length(lines2);
    
    lines1=addpointsnearby(lines1,pointlist1,sublinds1, pts1');
    lines2=addpointsnearby(lines2,pointlist2,sublinds2, pts2');

    simL=zeros(len1,len2);
    simR=zeros(len1,len2);

  disp(' Calculating similarities between line neighborhoods...');
    for i=1:len1
        t=lines1(i);
        for j=1:len2
          [simL(i,j),simR(i,j)]=distline(t,lines2(j));
        end 
    end

    k=[];
    for i=1:len1
        for j=1:len2
            if simL(i,j)>0.95&&( simL(i,j)==max(simL(i,:)) && simL (i,j) == max(simL(:,j)) )
                k=[k;[i,j]];
                break;
            end
        end
    end
    simside1=ones(1,size(k,1));
    
    for i=1:len1
        for j=1:len2
            if simR(i,j)>0.95&&( simR(i,j)==max(simR(i,:)) && simR (i,j) == max(simR(:,j)) )
                k=[k;[i,j]];
                break;
            end
        end
    end
simeside1=[ simside1 2*ones(1,size(k,1))];

len=size(k,1);
votecan=zeros(len2,len1);

 disp(' Matching lines ...');
for i=1:len

    [p1,p2]=getHpoints1L(lines1(k(i,1)),lines2(k(i,2)),simeside1(i));
      
    if ~isempty(p1) 
        [F1,~,~] =  estimateGeometricTransform(p1,p2,'projective');
        plines = projline(F1.T,lines1);
        [ind11,ind12]=getgoodpair(plines,lines2,3);
        
        plines = projline(inv(F1.T),lines2);
        [ind22,ind21]=getgoodpair(plines,lines1,3);
        
        if isempty(ind11)||isempty(ind22)
             continue;
        end

        [indfinal]=intersect([ind11;ind12]',[ind21;ind22]','rows');
        
        if ~isempty(indfinal)
            indfinal=indfinal';
            ind1=indfinal(1,:);
            ind2=indfinal(2,:);
        else
            ind1=[];ind2=[];
        end

         if simeside1(i)==1 
             v=simL(k(i,1),k(i,2));
         elseif  simeside1(i)==2 
              v=simR(k(i,1),k(i,2));
         end
 
        votecan(ind2+(ind1-1)*len2)=votecan(ind2+(ind1-1)*len2)+v;
    end
end

   [num,ind]=sort(votecan,'descend');
    num=num(1,:);
    ind=ind(1,:);
    votecan=votecan';
    [num2,ind2]=sort(votecan,'descend');
    num2=num2(1,:);
    ind2=ind2(1,:);
    k=[];
    for i=1:length(ind)
        if i==ind2(ind(i)) && (num(i) > 0.9 && num2(ind(i))> 0.9)
            k=[k,i];
        end
    end
%    nummatch=length(k);
%     draw(img1,lines1(k),1:nummatch,num2str(nummatch));   
%     draw(img2,lines2(ind(k)),1:nummatch,num2str(nummatch));  
    linestruct1 = lines1(k);
    linestruct2 = lines2(ind(k));
    linematch1 = zeros(length(linestruct1),4); linematch2 = linematch1;
    for i=1:length(linestruct1)
        linematch1(i,1:2)= linestruct1(i).point1;  % [x1, y1]
        linematch1(i,3:4)= linestruct1(i).point2;  % [x2, y2]
        linematch2(i,1:2)= linestruct2(i).point1;  % [x1, y1]
        linematch2(i,3:4)= linestruct2(i).point2;  % [x2, y2]
    end
    
    %---- delete outliers of line match
    linematch1 = refineLine(linematch1, img1);
    linematch2 = refineLine(linematch2, img2);
    [linematch1, linematch2] = linesDelete(linematch1, linematch2, pts1, pts2);
    
end

% function draw(I,lines,orders,name)
% 
% %I=imread(img);
% len=length(orders);
% 
%     figure, imshow(I),hold on 
%     if exist('name','var') 
%          title(name);
%     end
%         for k = 1:len
%             if orders(k)~=0
%                 xy = [lines(orders(k)).point1; lines(orders(k)).point2];
%                 plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red');
%                 text((xy(1,1)+xy(2,1))/2,(xy(1,2)+xy(2,2))/2,num2str(k));
%             end   
%         end
%     hold off;
% 
% end
