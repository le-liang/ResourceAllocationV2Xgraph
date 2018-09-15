%input w is the weight vector of the edges, and it is also a M*F*N matrix.
function [Z,fval,approValue]=solution3DMatching(w)
[M,F,N]=size(w);
A=zeros(M+F+N,M*F*N);

for a=1:1:M
    order=1;
    for i=1:1:M
        for j=1:1:F
            for k=1:1:N
                if i==a
                    A(a,order)=1;
                end
                order=order+1;
            end
        end
    end
end

for b=1:1:F
    order=1;
    for i=1:1:M
        for j=1:1:F
            for k=1:1:N
                if j==b
                    A(M+b,order)=1;
                end
                order=order+1;
            end
        end
    end
end

for c=1:1:N
    order=1;
    for i=1:1:M
        for j=1:1:F
            for k=1:1:N
                if k==c
                    A(M+F+c,order)=1;
                end
                order=order+1;
            end
        end
    end
end

b=ones(1,M+F+N);
f=-threeToOne(w);
Aeq=[];
Beq=[];
ub=[];
lb=zeros(1,M*F*N);

% options = optimoptions('linprog','Algorithm','simplex');
options = optimoptions('linprog','Algorithm','dual-simplex','display','off');
% options = optimoptions('linprog','Algorithm','interior-point');

x0=[];
[x,fval,exitflag,output] = linprog(f,A,b,Aeq,Beq,lb,ub,x0,options);
fval=-x'*f';
exitflag;
X=oneToThree(x,M,F,N);
Z=D3matching(X,w);

approValue=sum(sum(sum(Z.*w)));
end

function one=threeToOne(three)
[M,F,N]=size(three);
one=zeros(1,M*F*N);
order=1;
for i=1:1:M
    for j=1:1:F
        for k=1:1:N
            one(order)=three(i,j,k);
            order=order+1;
        end
    end
end
end


function three=oneToThree(one,M,F,N)
three=zeros(M,F,N);
order=1;
for i=1:1:M
    for j=1:1:F
        for k=1:1:N
            three(i,j,k)=one(order);
            order=order+1;
        end
    end
end
end


%input X is a M*F*N matrix, and it is a basic solution with each entry's
%value between 0 and 1. For example, X(m,f,n) denotes the value of the
%edge(m,f,n).

%input w is the weight vector of the edges, and it is also a M*F*N matrix.

%output Z is a M*F*N matrix, and all entries of it are 0 or 1.
function Z=D3matching(X, w)
dimension=3;
[m,f,n]=size(w);
Inside=ones(m,f,n);
F=zeros(m,f,n);

Sequence=zeros(1,3*m*f*n);

order=1;
while (order<=(m*f*n))
    Y=X.*Inside;
    
    for j=1:1:m
        for k=1:1:f
            for l=1:1:n
                if Inside(j,k,l)==1
                    sumValue=sumOfAllNeighbor(Y,j,k,l);
                    if sumValue<=(dimension-1)
                        F(j,k,l)=order;
                        
                        Sequence(3*order-2)=j;
                        Sequence(3*order-1)=k;
                        Sequence(3*order)=l;
                        
                        order=order+1;
                        Inside(j,k,l)=0;
                        Y(j,k,l)=0;
                    end
                end
            end
        end
    end
end
%Now, we get F and Sequence, which gives us the order of edges of hypergraph.
%preprocess F and w, such that F and w would both be zero, or not
judgeMatrix=(w~=zeros(m,f,n));
F=F.*judgeMatrix;
%Now start LocalRatio
Result=LocalRatio(F,w,Sequence);
Z=zeros(m,f,n);
[~,sizeResult]=size(Result);
yui1=zeros(1,m);
yui2=zeros(1,f);
yui3=zeros(1,n);
for i=3:3:sizeResult

    label1=Result(i-2);
    label2=Result(i-1);
    label3=Result(i);
    yui1(label1)=1;
    yui2(label2)=1;
    yui3(label3)=1;
    Z(label1,label2, label3)=1;
end
for i=1:1:m    
    for j=1:1:f        
        for k=1:1:n
            
         if yui1(i)==0
             if yui2(j)==0
        if yui3(k)==0    
            if w(i,j,k)>=0
               Z(i,j,k)=1; 
               yui1(i)=1;
               yui2(j)=1;
               yui3(k)=1;
            end           
        end
             end
         end
         
        end        
    end    
end

end

function Result=LocalRatio(F, w, Sequence)%F and w would both be zero, or not
%Result is the edges in matching
%Result=[3,3,3];

[wM,wF,wN]=size(w);
[~,si]=size(Sequence);
for i=3:3:si
    m=Sequence(i-2);
    f=Sequence(i-1);
    n=Sequence(i);
    if F(m,f,n)~=0
        value=w(m,f,n);
        Check=zeros(wM,wF,wN);
        
        %for m
        for a=1:1:wF
            for b=1:1:wN
                if Check(m,a,b)==1
                    
                else
                    Check(m,a,b)=1;
                    if w(m,a,b)>value
                        w(m,a,b)=w(m,a,b)-value;
                    else
                        w(m,a,b)=0;
                        F(m,a,b)=0;
                    end
                end
            end
        end
        
        % for f
        for a=1:1:wM
            for b=1:1:wN
                if Check(a,f,b)==1
                    
                else
                    Check(a,f,b)=1;
                    if w(a,f,b)>value
                        w(a,f,b)=w(a,f,b)-value;
                    else
                        w(a,f,b)=0;
                        F(a,f,b)=0;
                    end
                end
            end
        end
        
        
        %for n
        for a=1:1:wM
            for b=1:1:wF
                if Check(a,b,n)==1
                    
                else
                    Check(a,b,n)=1;
                    if w(a,b,n)>value
                        w(a,b,n)=w(a,b,n)-value;
                    else
                        w(a,b,n)=0;
                        F(a,b,n)=0;
                    end
                end
            end
        end
        
        %Now, we judge whether or not we need to continue to next LocalRatio
        Kmatrix=(F~=zeros(wM,wF,wN));
        sumKmatrix=sum(sum(sum(Kmatrix)));
        if sumKmatrix>0
            Sequence=Sequence(i+1:si);
            
            Result2=[m,f,n];
            Result1=LocalRatio(F,w,Sequence);
            
            %We need to judge whether or not Result2 and Result1 could form
            %a matching
            if Judge(Result1,Result2)==1
                [~,sizeRe]=size(Result1);
                Result=zeros(1,sizeRe+3);
                Result(1:3)=Result2;
                Result(4:sizeRe+3)=Result1;
            else
                Result=Result1;
            end
            
        else
            Result=[m,f,n];
        end
        
        %Now End
        break;
    end
end
end


%We need to calculate X(Neighbor(m,f,n)).
function sumValue=sumOfAllNeighbor(X, m, f, n)
sum1=sum(sum(X(m,:,:)));
sum2=sum(sum(X(:,f,:)));
sum3=sum(sum(X(:,:,n)));
sum4=sum(X(m,f,:));
sum5=sum(X(m,:,n));
sum6=sum(X(:,f,n));
sum7=X(m,f,n);
sumValue=sum1+sum2+sum3-sum4-sum5-sum6+sum7;
end

%We need to judge whether or not Result1 and Result2 could form a matching
function isMatching=Judge(Result1, Result2)
[~,long]=size(Result1);
isMatching=1;
for i=3:3:long
    if Result2(1)==Result1(i-2)
        isMatching=0;
        break;
    end
    if Result2(2)==Result1(i-1)
        isMatching=0;
        break;
    end
    if Result2(3)==Result1(i)
        isMatching=0;
        break;
    end
end
end