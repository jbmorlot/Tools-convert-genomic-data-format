function mat2bed(vect,filename,res)
disp(['Output file:' filename,'.bed'])
fileID = fopen([filename '.bed'],'w');

chrmpos = [0,249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566];
CSchrmpos = cumsum(chrmpos);

BCSchrmpos = cumsum(floor(chrmpos./res));
chrm=0;
for i=1:length(vect)
    
    if i>BCSchrmpos(chrm+1)
        chrm = chrm+1;
        if chrm>1
           toc, 
        end
        disp(['Chromosom ' num2str(chrm)])
        tic;
    end
    if floor(vect(i))>0
        pos = ((i-1)-BCSchrmpos(chrm))*res;

        %[chrm,startp,endp] = bin2epigenom(peaksidx(i),peaksidx(i),res);
        fprintf(fileID,'chr%d\t%d\t%d\t%d\n',chrm,pos,pos+res,floor(vect(i)));
    end
end

fclose(fileID);


end
