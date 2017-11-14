function [] = writeSolution(caseNumber,ue,u1,u2,u3,u0,c3,Ga3,G3,error0,error2,error3, L2ele0, L2ele2,L2ele3,H0norm,H2norm,H3norm)
casenumber = string(caseNumber);
filename1 = strcat('./udata/case',casenumber,'.m');
filename2 = strcat('./udata/case',casenumber,'dynamic.m');
filename3 = strcat('./udata/case',casenumber,'error.m');

fileID = fopen(filename1,'w');
string1 = strcat('function [u1,u2,u3,ue,u0] = case',casenumber,'() \n');
fprintf(fileID, string1);

fprintf(fileID, 'u1 = [');
fprintf(fileID,'%6.10f \n',u1);
fprintf(fileID, '];\n');

fprintf(fileID, 'u2 = [');
fprintf(fileID,'%6.10f \n',u2);
fprintf(fileID, '];\n');

fprintf(fileID, 'u3 = [');
fprintf(fileID,'%6.10f \n',u3);
fprintf(fileID, '];\n');

fprintf(fileID, 'ue = [');
fprintf(fileID,'%6.10f \n',ue);
fprintf(fileID, '];\n');

fprintf(fileID, 'u0 = [');
fprintf(fileID,'%6.10f \n',u0);
fprintf(fileID, '];\n');
fprintf(fileID, 'end \n');
fclose(fileID);

fileID = fopen(filename2,'w');
string2 = strcat('function [c1,Ga,GIter] = case',casenumber,'dynamic() \n');
fprintf(fileID, string2);
fprintf(fileID, 'c1 = [');
fprintf(fileID,'%6.10f \n',c3);
fprintf(fileID, '];\n');

fprintf(fileID, 'Ga = [');
fprintf(fileID,'%6.10f \n',Ga3);
fprintf(fileID, '];\n');

fprintf(fileID, 'GIter = [');
fprintf(fileID,'%6.10f \n',G3);
fprintf(fileID, '];\n');
fprintf(fileID, 'end \n');
fclose(fileID);

fileID = fopen(filename3,'w');
string3 = strcat('function [error0,error2,error3,H2norm,H3norm,H0norm,L2ele0,L2ele2,L2ele3] = case',casenumber,'error() \n');
fprintf(fileID, string3);
fprintf(fileID, 'error0 = ');
fprintf(fileID,'%6.16f;\n',error0);

fprintf(fileID, 'error2 = ');
fprintf(fileID,'%6.16f; \n',error2);

fprintf(fileID, 'error3 = ');
fprintf(fileID,'%6.16f; \n',error3);

fprintf(fileID, 'H2norm = ');
fprintf(fileID,'%6.16e;\n',H2norm);

fprintf(fileID, 'H3norm = ');
fprintf(fileID,'%6.16e; \n',H3norm);

fprintf(fileID, 'H0norm = ');
fprintf(fileID,'%6.16e; \n',H0norm);

fprintf(fileID, 'L2ele0 = [');
fprintf(fileID,'%6.16f \n',L2ele0);
fprintf(fileID, '];\n');

fprintf(fileID, 'L2ele2 = [');
fprintf(fileID,'%6.16f \n',L2ele2);
fprintf(fileID, '];\n');

fprintf(fileID, 'L2ele3 = [');
fprintf(fileID,'%6.16f \n',L2ele3);
fprintf(fileID, '];\n');
fprintf(fileID, 'end \n');
fclose(fileID);
end