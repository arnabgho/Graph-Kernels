function [result,runtime] = fastgraphkernelmatrix(graph, options) 
% calculate results for pairs of graphs


% preprocessing

% if any(strcmp('restrict',options))
%   random_indices=randperm(size(graph,2))
%   graph=graph(random_indices)
t=tic

if any(strcmp('edgelabels',options))
display('edgelabels');
%edgelabel2 = zeros(size(am2,1),size(am2,2));

for j = 1:size(graph,2)

  edgelabel1 = zeros(size(graph(j).am,1),size(graph(j).am,2));

for i = 1:size(graph(j).el.values,1)
if (graph(j).el.values(i,1) > 0) && (graph(j).el.values(i,2) > 0)
edgelabel1(graph(j).el.values(i,1),graph(j).el.values(i,2))= graph(j).el.values(i,3); 
edgelabel1(graph(j).el.values(i,2),graph(j).el.values(i,1))= graph(j).el.values(i,3);
end
end

am1x = full(graph(j).am .* edgelabel1);

%i1 = intersect(graph(j).el.values(:,3),graph(j).el.values(:,3))';
graph(j).am1 = filteram(am1x);%,i1);
%i2 = intersect(g2.el.values(:,3),g2.el.values(:,3))';
%am2 = filteram(am2x,i2);
end
end

graph2 = graph;



  for i = 1:size(graph,2)
    % i
    for j = 1:size(graph2,2)
      %j
%try
      result(i,j) =  fastkernel(graph(i), graph2(j),options);
%catch
%result(i,j) = Inf;
%end
      end
  end
  % runtime=cputime-t   
  runtime= toc(t)
% end    
function result = fastkernel(g1, g2,options)
% General variables
result = 0;
%lambda = 0.0001; %MUTAG, ENZYMES, PTC
lambda = 1e-6;

am1 = g1.am;
am2 = g2.am; 

if (any(strcmp('fix',options)) || any(strcmp('dir',options)) || ...
    any(strcmp('cg',options)) ||  any(strcmp('syl',options)) ||  any(strcmp('linear',options)) ||  any(strcmp('sparse',options)) ||  any(strcmp('fixnovec',options)) )
%display('multiplied')
%WW = kron(am1,am2);
am2 = am2 * lambda;
end

% #################################################
% choose edge kernel
% #################################################

% for linear kernel
if any(strcmp('linear',options))
%display('linear');

%measuretime = cputime;
edgelength1 = zeros(size(am1,1),size(am1,2));
edgelength2 = zeros(size(am2,1),size(am2,2));

for i = 1:size(g1.el.values,1)
if (g1.el.values(i,1) > 0) && (g1.el.values(i,2) > 0)
edgelength1(g1.el.values(i,1),g1.el.values(i,2))= g1.el.values(i,10);
edgelength1(g1.el.values(i,2),g1.el.values(i,1))= g1.el.values(i,10);
end
end

for i = 1:size(g2.el.values,1)
if (g2.el.values(i,1) > 0) && (g2.el.values(i,2) > 0)
edgelength2(g2.el.values(i,1),g2.el.values(i,2))= g2.el.values(i,10);
edgelength2(g2.el.values(i,2),g2.el.values(i,1))= g2.el.values(i,10);
end
end

am1 = am1 .* edgelength1;
am2 = am2 .* edgelength2;

%runtime_preprocessing = cputime - measuretime

end

% for delta kernel

if any(strcmp('edgelabels',options))
%display('edgelabels');
edgelabel1 = zeros(size(am1,1),size(am1,2));
edgelabel2 = zeros(size(am2,1),size(am2,2));

for i = 1:size(g1.el.values,1)
if (g1.el.values(i,1) > 0) && (g1.el.values(i,2) > 0)
edgelabel1(g1.el.values(i,1),g1.el.values(i,2))= g1.el.values(i,3); 
edgelabel1(g1.el.values(i,2),g1.el.values(i,1))= g1.el.values(i,3);
end
end

for i = 1:size(g2.el.values,1)
if (g2.el.values(i,1) > 0) && (g2.el.values(i,2) > 0)
edgelabel2(g2.el.values(i,1),g2.el.values(i,2))= g2.el.values(i,3);
edgelabel2(g2.el.values(i,2),g2.el.values(i,1))= g2.el.values(i,3);
end
end

am1x = full(am1 .* edgelabel1);
am2x = full(am2 .* edgelabel2);

i1 = intersect(g1.el.values(:,3),g1.el.values(:,3))';

am1 = filteram(am1x,i1);

i2 = intersect(g2.el.values(:,3),g2.el.values(:,3))';

am2 = filteram(am2x,i2);

%for i = 1:size(am2,2)
%am2(i).am =  am2(i).am;
%end
end


% #################################################
% choose graph kernel computation technique
% #################################################

if any(strcmp('fix',options))
%display('fix');
x = afixiterate(am1,am2,ones(size(am1,1)*size(am2,1),1),1e-6,20);
result = sum(sum(x));
end

if any(strcmp('fixnovec',options))
%display('fix_novec');
W = sparse(kron(am2',am1));
x = novecfixiterate(W,ones(size(am1,1)*size(am2,2),1),1e-6,20);
result = sum(sum(x));
end

if any(strcmp('cg',options))
%display('cg');
[x,rubbish] = pcg(@(x)smt(x,am1,am2),ones(size(am1,1)*size(am2,1),1),1e-6,20);
result =  sum(sum(x)); 
end
 
if any(strcmp('syl',options))
%display('syl');
am1 = full(am1);
am2 = full(am2);
X =dlyap(am2,am1,ones(size(am2,1),size(am1,1)));
result = sum(sum(X));
end

%%%%

if any(strcmp('dir',options))
display('dir');
W = kron(am1,am2);

%try
zresult = inv(speye(size(W)) - W) ;
%catch
%display('Out of memory');
%zresult = 0;
%end

result = ones(1,size(W,1)) * zresult * ones(1,size(W,1))';
end

%%%%

if any(strcmp('sparse',options))
%display('sparse');
W = sparse(kron(am1,am2));
[Wsize1, Wsize2] = size(W);

zresult = sparse((speye(Wsize1,Wsize2) - W)) \ ones(Wsize1,1);
result = ones(1,Wsize1) * zresult;
clear zresult W
end

%###########################################################################
%
%##########################################################################

% #################################################
% choose graph kernel computation technique for delta
% #################################################


%%% new vec trick
if any(strcmp('fixdelta',options))
%display('fixdelta');

am1 = g1.am1;
am2 = g2.am1;



intersection = intersect(g1.el.values(:,3),g2.el.values(:,3))';

x = afixiteratefilter(am1,am2,ones(size(g1.am,1)*size(g2.am,1),1),1e-6,20,lambda,intersection);
result = sum(sum(x));
end


%%% new vec trick
if any(strcmp('cgdelta',options))
%display('cgdelta');
am1 = g1.am1;
am2 = g2.am1;
% determine  x from (I-W)*x = x_0
[x,rubbish] = pcg(@(x)smtfilter(x,am1,am2,lambda),ones(size(g1.am,1)*size(g2.am,1),1),1e-6,20);
% x_0' * x
result =  sum(sum(x));
end




%%% KPA 
if any(strcmp('syldelta',options))
display('syldelta');
am1 = g1.am1;
am2 = g2.am1;

W=zeros(size(g1.nl.values,1) * size(g2.nl.values,1),size(g1.nl.values,1) * size(g2.nl.values,1));
for i = 1:4
W = W + kron(am2(i).am',am1(i).am); 
end
W = W * lambda;

[T,S] = kpa(full(W),size(g1.am,1),size(g1.am,2));
X=dlyap(S,T',ones(size(S,1),size(T,1)));

%X =dlyap(am2,am1,ones(size(am2,1),size(am1,1)));
result = sum(sum(X));
end



%%% direct 
if any(strcmp('dirdelta',options))
% display('dirdelta')
am1 = g1.am;
am2 = g2.am;
W=sparse(zeros(size(g1.nl.values,1) * size(g2.nl.values,1),size(g1.nl.values,1) * size(g2.nl.values,1)));

for i = 1:4
% W = W + kron(am1(i).am,am2(i).am); 
  W = W + kron(am1,am2); 
end
%[V,D]=eig(W);
%largesteigenvalue = max(diag(D));
%sum(sum(W));
%W
%nnz(W);
zresult = inv(eye(size(W,1)) - lambda*W) ;
result = ones(1,size(W,1)) * zresult * ones(1,size(W,1))';
end



%%% Shortest path kernel
if any(strcmp('sparsedelta',options))
%display('sparsedelta');
am1 = g1.am;
am2 = g2.am;
W=sparse(zeros(size(am1,1) * size(am2,1),size(am1,1) * size(am2,1)));


for i = 1:4
additive_term =   kron(am1,am2x);
if (isempty(additive_term) == 0) 
W = sparse(W + additive_term);
end
end

[Wsize1, Wsize2] = size(W);
zresult = sparse((speye(Wsize1,Wsize2) - (W*lambda)) \ ones(Wsize1,1));
result = ones(1,Wsize1) * zresult;
clear zresult W
end


%%% shortest path kernel

if any(strcmp('sp',options))
% display('sp');
sp1 = floydwarshall(am1);
sp1(find(sp1 == Inf)) = 0;

sp2 = floydwarshall(am2);
sp2(find(sp2 == Inf)) = 0;

maxlen = 0 + min(max(max(sp1)),max(max(sp2)));
result = 0;
for i = 1:maxlen
result = result + sum(sum(sp1==i))*sum(sum(sp2==i));
end
end
 


%%% shortest path kernel

if any(strcmp('productgraph',options))
%display('pg');

[n1, nothing]= size(g1.nl.values);
[n2, nothing]= size(g2.nl.values);


g1.v(:,2) = repmat(g1.nl.values,n2,1);
g1.v(:,1) = repmat([1:n1]',n2,1);

g2.v(:,2) = vec(repmat(g2.nl.values',n1,1));
g2.v(:,1) = vec(repmat([1:n2],n1,1));

% g2_temp2=repmat(g2.nl.values',n1,1);
% g2_temp1=repmat([1:n2],n1,1);

% g2.v(:,2)=reshape(g2_temp2.',1,[]);
% g2.v(:,1)=reshape(g2_temp1.',1,[]);

matches = [g1.v(:,1), g2.v(:,1), (g1.v(:,2) -g2.v(:,2))];

[same, notused] = find(matches(:,3)==0);

W = sparse(g1.am(g1.v(same,1),g1.v(same,1)) .* g2.am(g2.v(same,1), ...
						  g2.v(same,1))) * lambda;

[Wsize1, Wsize2] = size(W);

zresult = sparse((speye(Wsize1,Wsize2) - W) \ ones(Wsize1,1));
result = ones(1,Wsize1) * zresult;

end


% floydwarshall: function description
function [sp] = floydwarshall(am)
  sp = full(am);
  sp(find(sp == 0)) = Inf;
  for k = 1:size(sp,1)
    for i = 1:size(sp,1)
      for j = 1:size(sp,1)
        sp(i,j)=min(sp(i,j),sp(i,k)+sp(k,j));
      end
    end  
  end
% end

function v = vec(M)
[m,n] = size(M);
v = reshape(M,m*n,1);
% end

function v = invvec(M,m,n)
v = reshape(M,m,n);

% end