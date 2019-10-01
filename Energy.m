function [E,amino,new]=Energy(e,x,y)

%{
Function to calculate the energy of the native state of the protein.
Arguments:
e - Energy of a single non-covalent interaction
x and y - coordinates of the native state protein
returns:
E - The calculated total energy of the native state
new - A cell array which consists of the coordinates
amino - A 16x3 cell array which keeps track of the native state
interactions
%}

l=length(x);
c=0;
new=cell(16,1);
amino=cell(16,3);
for i = 1:l
    new{i}=[x(i) y(i)];
end
for i = 1:l
    q=[x(i)+1 y(i)];
     isQ = cellfun(@(x)isequal(x,q),new);
     is1 = find(isQ);
    r=[x(i) y(i)+1];
     isR = cellfun(@(x)isequal(x,r),new);
     is2 = find(isR);
    s=[x(i)-1 y(i)];
     isS = cellfun(@(x)isequal(x,s),new);
     is3 = find(isS);
    t=[x(i) y(i)-1];
     isT = cellfun(@(x)isequal(x,t),new);
     is4 = find(isT);
    k=1;

        if ~isempty(is1) && is1~=i-1 && is1~=i+1
            c=c+1;
            amino{i,k}=find(isQ==1);        %store the index of the residue it is interacting with
            k=k+1;
        end
       if ~isempty(is2) && is2~=i-1 && is2~=i+1
            c=c+1;
            amino{i,k}=find(isR==1);
            k=k+1;
       end
      if ~isempty(is3) && is3~=i-1 && is3~=i+1
            c=c+1;
            amino{i,k}=find(isS==1);
            k=k+1;
      end
      if ~isempty(is4) && is4~=i-1 && is4~=i+1
            c=c+1;
            amino{i,k}=find(isT==1);
      end
        
      
end
            
E=e*c/2;    %to avoid double counting
end

    