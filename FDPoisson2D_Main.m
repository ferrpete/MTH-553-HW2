clear all

ax=0;
bx=11;
ay=0;
by=28;

Mx=176;
My=448;
hx=bx/Mx;
hy=by/My;

x=hx*(0:1:Mx);
y=hy*(0:1:My);
[xx,yy]=meshgrid(x,y); %% ’’grid’’ has Mx columns and My rows
%%
myin = zeros(size(xx'));
for i=1:Mx+1
    for j=1:My+1
        myin(i,j) = myind(i,j,Mx,My);
    end
end
% %% show the numbering provided by index
% zz,
% %% show the numbering as intuition suggests
% flipud(zz'),
% %%
% surfc(xx,yy,zz');title('Index of grid points');pause;
% %% show the interpretation of grid
% surfc(xx,yy,xx/hx+yy/hy);title('Affine function of grid-points');
% %% draw just the mesh
% mesh(xx,yy,zz');

Nall = (Mx+1)*(My+1);
F = sparse(Nall,1);
A = sparse(Nall,Nall);

Uexact=0*F;

interior_nodes = [];
exterior_nodes = [];
load('SEdge.txt')
xEdge = SEdge(:,1);
yEdge = SEdge(:,2);

for i=1:Mx+1
    for j=1:My+1
        me = myin(i,j);
        [in, on] = inpolygon(x(i),y(j),xEdge,yEdge);
        if in == 1 && on == 0
            interior_nodes = [interior_nodes,me];
            F(me)= rhsfun(x(i),y(j));
            Uexact(me) = exactfun(x(i),y(j));
        else
            exterior_nodes = [exterior_nodes,me];
            F(me)= 0;
        end
    end
end

%%%
for i=2:Mx
    for j=2:My
        me = myin(i,j);
        mel=myin(i-1,j);
        mer=myin(i+1,j);
        meb=myin(i,j-1);
        met=myin(i,j+1);
        %%
        if ismember(me, interior_nodes)
            A(me,me)=2/hx^2+2/hy^2;
        end
        %% neighbors might be on the boundary.
        A(me,mel)=-1/hx^2;
        
        if ismember(mel,interior_nodes)
            A(mel,me)=-1/hx^2;
        end
        
        A(me,mer)=-1/hx^2;
        
        if ismember(mer,interior_nodes)
            A(mer,me)=-1/hx^2;
        end
        
        A(me,meb)=-1/hy^2;
        
        if ismember(meb,interior_nodes)
            A(meb,me)=-1/hy^2;
        end
        
        A(me,met)=-1/hy^2;
        
        if ismember(met,interior_nodes)
            A(met,me)=-1/hy^2;
        end
    end
end
%% boundary conditions imposed explicitly

for i=1:Mx+1
    for j=1:My+1
        me = myin(i,j);
        if ismember(me, exterior_nodes)
            A(me,me) = 1;
        end
    end
end
% i=1;
% for j=1:My+1
%     me = myin(i,j);
%     A(me,me)=1;
% end
% 
% i = Mx+1;
% for j=1:My+1
%     me = myin(i,j);
%     A(me,me)=1;
% end
% 
% j = 1;
% for i=1:Mx+1
%     me = myin(i,j);
%     A(me,me)=1;
% end
% 
% j=My+1;
% for i=1:Mx+1
%     me = myin(i,j);
%     A(me,me)=1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = A\F;
%% unpack the solution
for i=1:Mx+1
    for j=1:My+1
        me=myin(i,j);
        if ismember(me,exterior_nodes)
            Uplot(i,j) = 0;
        else
            Uplot(i,j)=U(me);
        end
        Uexactplot(i,j)=Uexact(me);
        
    end
end
%% plot and check convergence
surf(xx,yy,Uplot');title('Numerical solution');
fprintf('Error =%g\n',norm(Uexact-U,inf));


%%%%%%%%%%%%%%%%%%%%%%%
function ind=myind(i,j,Mx,My)

ind=(j-1)*(Mx+1) + i;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%