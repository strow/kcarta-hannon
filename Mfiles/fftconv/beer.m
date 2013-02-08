function apod=beer(N,MD);
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3
%

beer=zeros(N,1);
beer(1)=1;
for i=2:N/2
        if i <= MD
                beer(i)=(1-((i-1)/MD)^2)^2;
        else
                beer(i)=0;
        end
end
beer(N/2+1:N)=flipud(beer(1:N/2));

apod=beer;
