function fI=bind_fxn(h,i,j)
%Sam Furness
%2/18/19
%Takes concentration values for mRNA or protein and 
    %calculates hill function. I is the effector and J is the target.
KI1=0.3;%mM
nI1=1.5;
K12=1000;%nmol/gDW
n12=1.5;
K13=1000;%nmol/gDW
n13=1.5;
K23=10000;%nmol/gDW
n23=10;
if i==1 %effector
    if j==1 %target
        fI=(h^nI1)/(KI1^nI1+h^nI1);
    elseif j==2 
        fI=(h^n12)/(K12^n12+h^n12);
    else
        fI=(h^n13)/(K13^n13+h^n13);
    end
else
    fI=(h^n23)/(K23^n23+h^n23); 
end
    




end