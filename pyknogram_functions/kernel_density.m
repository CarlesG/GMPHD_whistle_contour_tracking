function  [P,id]  = kernel_density( FW, ndraw, BW, BWoverlap )
%   
%   This function obtains the possible whistle points from the pyknogram
%   based on a point cloud filterning using the Gausian kernel K(x)
%   referred as the Parzen window technique (Parzen, 1962).
%   Inspiered in the work: 
%       Kernel Density Filtering for Noisy Point Clouds in One Step,
%       M.A. Brophy, S.S. Beauchemin, J.L. Barron, IMVIP 2015
%
%   [P,id]  = kernel_density( FW, ndraw, BW, BWoverlap )
%

%radius=(BW*0.5*(100-BWoverlap)/100);
h=(BW*0.5*(100-BWoverlap)/100);

density=zeros(size(FW));

for k=1:size(FW,1)
    for l=1:size(FW,2)
        di=(FW(k,l)-FW(k,:))/h;
        K=abs(di)<1;  % Rect Kernel
        fxi=1/(h*length(di))*sum(K);
        density(k,l)=fxi;
    end
end

[ni,nj]=find(density>1/(h*length(di)));
id=sub2ind(size(FW),ni,nj);

P(1:length(ni))=struct('time',[],'freq',[],'ampl',[],'label',[],'done',[]); 

kk=num2cell(ndraw(ni));[P.time]=kk{:};
kk=num2cell(FW(id));[P.freq]=kk{:};

% figure(10); sh=scatter([P.time],[P.freq],'filled');
% sh.SizeData=15;


% density=zeros(size(FW));
% 
% d=2;  % Dimension
% h
% N=
% for i1=1:prod(size(FW))
%     [k,l]=ind2sub(size(FW),i1);
%     xi=[ndraw(k) FW(k,l)];
%     for i2=1:prod(size(FW))
%         [k,l]=ind2sub(size(FW),i2);
%         xj=[ndraw(k) FW(k,l)];
%         junk = (xi - xj);
%         euc_dist_squared=(junk * junk');
%     end
% end

end

