function inout=standVectorUp(inout)
   [y,~]=size(inout);
   if y==1
       inout=inout';
   end 
end