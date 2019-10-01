function E=EnergyNew(e,amino,x,y)

%{
Function to calculate the energy of the current conformation of the protein.
Arguments:
e - Energy of a single non-covalent interaction
x and y - coordinates of the native state protein
amino - A 16x3 cell array which keeps track of the native state
interactions 
Returns:
E - The calculated total energy of the current conformation
%}

l=length(x);
c=0;
new=cell(16,1);
for i = 1:l
    new{i}=[x(i) y(i)];       %Storing the positions of the residues in a cell array
end
for i = 1:l
    interactions=amino(i,:);
    orig=new{i,1};
    z1=interactions{1,1};
    z2=interactions{1,2};
    z3=interactions{1,3};
    if ~isempty(z1)             %checking if it is a native state interaction
        inter=new{z1,1};
        if norm((orig-inter),2)==1
            c=c+1;
        end
    end
    if ~isempty(z2)
        inter=new{z2,1};
        if norm((orig-inter),2)==1
            c=c+1;
        end
    end
    if ~isempty(z3)
        inter=new{z3,1};
        if norm((orig-inter),2)==1
            c=c+1;
        end
    end

end

E=e*c/2;       %To avoid double counting

end


    