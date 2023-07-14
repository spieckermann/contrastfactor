function Compliance = tensor2voigt(TCompliance)

% function tensor2voigt
%
% DESCRIPTION
% Transforms the compliance from the  the tensor (rank 4) notation to voigt matrix notation (6x6)

Compliance = zeros(6,6);

Compliance(1,1) = TCompliance(1,1,1,1);
Compliance(1,2) = TCompliance(1,1,2,2);
Compliance(1,3) = TCompliance(1,1,3,3);
Compliance(1,4) = TCompliance(1,1,2,3);
Compliance(1,5) = TCompliance(1,1,3,1);
Compliance(1,6) = TCompliance(1,1,1,2);
  
Compliance(2,1) = TCompliance(2,2,1,1);
Compliance(2,2) = TCompliance(2,2,2,2);
Compliance(2,3) = TCompliance(2,2,3,3);
Compliance(2,4) = TCompliance(2,2,2,3);
Compliance(2,5) = TCompliance(2,2,3,1);
Compliance(2,6) = TCompliance(2,2,1,2);
  
Compliance(3,1) = TCompliance(3,3,1,1);
Compliance(3,2) = TCompliance(3,3,2,2);
Compliance(3,3) = TCompliance(3,3,3,3);
Compliance(3,4) = TCompliance(3,3,2,3);
Compliance(3,5) = TCompliance(3,3,3,1);
Compliance(3,6) = TCompliance(3,3,1,2);
  
Compliance(4,1) = TCompliance(2,3,1,1);
Compliance(4,2) = TCompliance(2,3,2,2);
Compliance(4,3) = TCompliance(2,3,3,3);
Compliance(4,4) = TCompliance(2,3,2,3);
Compliance(4,5) = TCompliance(2,3,3,1);
Compliance(4,6) = TCompliance(2,3,1,2);
  
Compliance(5,1) = TCompliance(3,1,1,1);
Compliance(5,2) = TCompliance(3,1,2,2);
Compliance(5,3) = TCompliance(3,1,3,3);
Compliance(5,4) = TCompliance(3,1,2,3);
Compliance(5,5) = TCompliance(3,1,3,1);
Compliance(5,6) = TCompliance(3,1,1,2);
  
Compliance(6,1) = TCompliance(1,2,1,1);
Compliance(6,2) = TCompliance(1,2,2,2);
Compliance(6,3) = TCompliance(1,2,3,3);
Compliance(6,4) = TCompliance(1,2,2,3);
Compliance(6,5) = TCompliance(1,2,3,1);
Compliance(6,6) = TCompliance(1,2,1,2);

Compliance = Compliance;
endfunction

