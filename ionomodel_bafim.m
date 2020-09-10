function [altitude,ne,te,ti,coll,cO,cM2,cH]=ionomodel_bafim(heights,modinfo)
%
% This function only calls ionomodel_iri. The actual Baysian
% filtering is implemented in apriorimodel_bafim.m
%
% Ilkka Virtanen, 2019-2020
%

[altitude,ne,te,ti,coll,cO,cM2,cH] = ionomodel_iri(heights,modinfo);



end

