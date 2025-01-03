function [collt,collv,colllu,collru,coll,tspt,nort,norv,norlu,norru,nor]=Thorax_twolung_ventricle(NB1,NBH)
%%
[collt,tsptt,nort]=thorax(NB1);
[collv,tsptv,norv]=ventricle(NBH);
[colllu,tsptlu,norlu]=lung_left(NBH);
[collru,tsptru,norru]=lung_right(NBH);

coll=[collt;colllu;collru;collv];
tspt=[tsptt;tsptv;tsptlu;tsptru];
nor=[nort;norv;norlu;norru];

