function A = makehatch(hatch)

%  By Ben Hinkle, bhinkle@mathworks.com

%  This code is in the public domain. 


n = 6;
A=zeros(n);
switch (hatch)
 case '/'
  A = fliplr(eye(n));
 case '\'
  A = eye(n);
 case '|'
  A(:,1) = 1;
 case '-'
  A(1,:) = 1;
 case '+'
  A(:,1) = 1;
  A(1,:) = 1;
 case 'x'
  A = eye(n) | fliplr(diag(ones(n-1,1),-1));
 case '.'
  A(1:2,1:2)=1;
 otherwise
  error(['Undefined hatch pattern "' hatch '".']);
end
