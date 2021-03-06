function RR=Statistic_OneWayAnova(R,expGroups)
numgroups=length(expGroups);

% prepare data
RS=[];GroupIndex=[];numsam=0;
for g=1:numgroups
    a=R{g,1};
    if size(a,1)>1
        R{g,1}=a';
    end
    RS=[RS,R{g,1}];
    for sam=1:length(R{g,1})
        numsam=numsam+1;
        GroupIndex{1,numsam}=expGroups{g};
    end
end

% anova1
[p,tb1,stats]=anova1(RS,GroupIndex,'off');
Fstat = tb1{2,5};
DF=[tb1{2,3} tb1{3,3}];

RR.P=p;
RR.Fstat=Fstat;
RR.DF=DF;

% multicompare
rs=multcompare(stats,'display','off');
%RR.mc=rs;

pcg=[];m=0;
for g=1:numgroups-1
    
    for cg=g+1:numgroups
        m=m+1;
        a=find(rs(:,1)==g);
        b=find(rs(:,2)==cg);
        c = intersect(a,b);
        
        pcg(m)=rs(c,6);
    end
end

RR.cg=pcg;


    

            