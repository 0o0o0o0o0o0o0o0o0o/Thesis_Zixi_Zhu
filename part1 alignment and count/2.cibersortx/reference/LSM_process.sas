libname data xlsx 'D:\Thesis\reference dataset\immunoStates\41467_2018_7242_MOESM5_ESM.xlsx';

proc contents data=data._all_ nodetails; 
run;


proc print data=DATA.'IMMUNSTATES MATRIX'n (obs=5);
run; 

proc import datafile="D:\Thesis\reference dataset\immunoStates\41467_2018_7242_MOESM5_ESM.xlsx"
    out=result dbms=xlsx replace;
    sheet="IMMUNSTATES MATRIX";
    range="A2:U1000";  /* Adjust as needed */
    getnames=yes;
run;

proc print data=result (obs=5);
run; 


proc export data=result outfile="D:\Thesis\reference dataset\immunoStates\ISM.txt" dbms=dlm 
    replace;
    delimiter='09'x;
    putnames=no;
run;

