function I_element=findelement(element)
Q=load('GRI3.0\CH4.mat');
elements=Q.elements;
I_element=find(strcmpi(element,elements));

end