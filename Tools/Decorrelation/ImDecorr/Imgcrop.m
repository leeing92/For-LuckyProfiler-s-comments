function Imgcrop(Im,subsize)
    Height = size(Im,1);
    Width = size(Im,2);
    HeightNum = floor(Height/subsize) +1;
    WidthNum = floor(Width/subsize)+1;
    oI = zeros(subsize,subsize);
    
    n=0;
    for i=1:HeightNum
        for j= 1:WidthNum
            if i<HeightNum && j<WidthNum
                oI = Im((i-1)*subsize+1:i*subsize,(j-1)*subsize+1:j*subsize);
                n=n+1;
                mysavetif(oI,n,subsize);
            elseif i==HeightNum && j<WidthNum
                temp =Im((i-1)*subsize+1:Height,(j-1)*subsize+1:j*subsize);
                if ~isempty(temp)
                oI(1:size(temp,1),1:size(temp,2)) =temp;
                n=n+1;
                mysavetif(oI,n,subsize);
                end
            elseif i<HeightNum && j==WidthNum
                temp =Im((i-1)*subsize+1:i*subsize,(j-1)*subsize+1:Width);
                if ~isempty(temp)
                  oI(1:size(temp,1),1:size(temp,2)) =temp;
                  n=n+1;
                  mysavetif(oI,n,subsize);
                end
            elseif i==HeightNum && j==WidthNum
                temp =Im((i-1)*subsize+1:Height,(j-1)*subsize+1:Width);
                if ~isempty(temp)
                oI(1:size(temp,1),1:size(temp,2)) =temp;
                n=n+1;
                mysavetif(oI,n,subsize);
                end
            end
        end
    end
    

end

function mysavetif(oI,n,subsize)
   str = ['E:\Mine\LAB\Resolution\Resolution research\Fig2\ImDecorrBlock\crop1\subsize ',num2str(subsize),'_',num2str(n),'.tif'];
   imwrite(uint16(oI),str);
end