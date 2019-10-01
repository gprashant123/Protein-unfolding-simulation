function ProteinFold(e,kT,num_iter)

%{
Main function - ProteinFold()
This function initializes the native state of the protein and finds the energy of the current conformation, having arguments 
e (Energy of a single non-covalent interaction) and kT(product of k and temperature). Further, it simulates the unfolding of proteins using
Monte Carlo Method. A plot is displayed after each iteration to show the
current state of the protein. This function makes use of two functions
Energy() to calculate the energy of the native state protein and EnergyNew() to calculate the energy of the current state 
%}

%initializing native state of protein
% x=[6 5 4 3 3 4 5 6 6 5 5 4 3 3 4 5];
% y=[7 7 7 7 6 6 6 6 5 5 4 4 4 3 3 3];

x = [1 0 0 1 1 0 0 1 2 2 2 2 3 3 3 3];
y = [3 3 2 2 1 1 0 0 0 1 2 3 3 2 1 0];


EVals=zeros(1,num_iter);

for i=1:num_iter
    count=0;
    if i==1
       [Ei,amino,new]=Energy(e,x,y);
    else
        
    pos=randi(16,1);   %randomly choosing which residue to move
    if pos==1 || pos==16    %The N-terminus and C-terminus residue can only undergo corner move
        if pos==1
            coord=new{pos+1,1};
        else
            coord=new{pos-1,1};
        end
        q=[coord(1)+1 coord(2)];
        g1=q-coord;
        g2=new{pos,1}-coord;
        r=[coord(1) coord(2)+1];
        h1=r-coord;
        h2=new{pos,1}-coord;
        s=[coord(1)-1 coord(2)];
        i1=s-coord;
        i2=new{pos,1}-coord;
        t=[coord(1) coord(2)-1];
        j1=t-coord;
        j2=new{pos,1}-coord;
        isQ = cellfun(@(x)isequal(x,q),new);
        isR = cellfun(@(x)isequal(x,r),new);
        isS = cellfun(@(x)isequal(x,s),new);
        isT = cellfun(@(x)isequal(x,t),new);
        tmp1=x(pos);
        tmp2=y(pos);
        if isempty(find(isQ,1))==1 && g1*g2'==0
            x(pos)=q(1);
            y(pos)=q(2);
        elseif isempty(find(isR,1))==1 && h1*h2'==0
            x(pos)=r(1);
            y(pos)=r(2);
        elseif isempty(find(isS,1))==1 && i1*i2'==0
            x(pos)=s(1);
            y(pos)=s(2);
        elseif isempty(find(isT,1))==1 && j1*j2'==0
            x(pos)=t(1);
            y(pos)=t(2);
        end
        Ef=EnergyNew(e,amino,x,y);
        new{pos}=[x(pos) y(pos)];
        w=exp(-(Ef-Ei)/kT);                %w is the probability of making the move
        count=1;
        if w<1
            var=rand();
            if var>w                           %if the random variable is greater the w, do not make the move
                x(pos)=tmp1;
                y(pos)=tmp2;
                new{pos}=[x(pos) y(pos)];
                count=0;
            end
        end
             
        
    else                           %The residues in between can undergo crank shaft move
        tmp1=x(pos);
        tmp2=y(pos);
        coord1=new{pos-1,1};
        coord2=new{pos,1};
        coord3=new{pos+1,1};
        n1=coord1-coord2;
        n2=coord3-coord2;
        if n1*n2'==0                %checking if the j-1,j,j+1 residue are at right angle to each other
            z=n1+n2;
            xn=z(1)+coord2(1);
            yn=z(2)+coord2(2);  
            cond = cellfun(@(x)isequal(x,[xn yn]),new);
            if isempty(find(cond,1))==1                 %checking if the new position is not already occupied by a residue
                x(pos)=xn;
                y(pos)=yn;
                new{pos}=[x(pos) y(pos)];
                Ef=EnergyNew(e,amino,x,y);
                w=exp(-(Ef-Ei)/kT);
                count=1;
                if w<1
                    var=rand();
                    if var>w
                        x(pos)=tmp1;
                        y(pos)=tmp2;
                        new{pos}=[x(pos) y(pos)];                
                        count=0;
                    end
                end
            end
        end    
    end
    end
    
    if count==0
     Ef=Ei;
    else
    Ei=Ef;
    end
EVals(i)=Ef;    

%Visualization

plot(x,y,'--o','MarkerFaceColor',[0.5 0 1]) 
hold on
plot(x(1),y(1),'s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor',[0 1 0])   %N terminus
plot(x(16),y(16),'s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor',[1 0 0])    %C terminus
xlim([-6 10])
ylim([-6 9]) 
txt = ['Energy:' num2str(Ef) ' units'];
txt1= ['Iteration:' num2str(i)];  
text(-5,9.5,txt1);
text(7.5,9.5,txt);
grid on               %for showing grid lines
grid minor
pause(0.07) %pausing the screen
hold off


% disp(i)
end

% plot(EVals,'-');
% xlabel('No. of iterations')
% ylabel('Energy of conformation')

