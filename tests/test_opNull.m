function test_suite = test_opNull
%test_opNull  Unit tests for the opNull operator
initTestSuite;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seed = setup
   seed = rng('default');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function test_opNull_noOperations(seed)
   % opNull must not result in changed array
   n = randi(100);
   x = randn(n,1) + 1i*randn(n,1);
   A = opNull(n);
   
   assertEqual( A*x, x )
   assertEqual( A'*x, x )
   assertEqual( A\x, x )
   assertEqual( A'\x, x )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function test_opNull_noWarning(seed)
   % opNull does not throw warning when applied to unknown types
   n = randi(100);
   x_Char = repmat('a',[n,1]);
   A = opNull(n);
   
   % ensure mtimes does not produce a warning
   % turn off warning verbosity for this stage of the test
   warning_settings = warning('off','all');
   
   % (first erase warning history with a fake warning)
   warning('opSpot:test:fakeWarning','fake warning');
   % now do a mtimes and check last warning
   y = A*x_Char;
   [lastWarning, lastWarningID] = lastwarn();
   
   % resume warning verbosity
   warning(warning_settings)
   
   assertTrue(not(strcmp(lastWarningID,'opSpot:mtimes:unsupportedType')))
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function test_opNull_noTypeChange(seed)
   % opNull result must be same as input type
   n = randi(100);
   x_cplxD = randn(n,1,'double') + 1i*randn(n,1,'double');
   x_cplxS = randn(n,1,'single') + 1i*randn(n,1,'single');
   x_D = randn(n,1,'double');
   x_S = randn(n,1,'single');
   x_I8 = randi(1000,n,1,'int8');
   x_Char = repmat('a',[n,1]);
   
   A = opNull(n);
   
   assertTrue(isa(A * x_cplxD, 'double')); assertTrue(not(isreal(A * x_cplxD)));
   assertTrue(isa(A * x_cplxS, 'single')); assertTrue(not(isreal(A * x_cplxS)));
   assertTrue(isa(A * x_D, 'double'));
   assertTrue(isa(A * x_S, 'single'));
   assertTrue(isa(A * x_I8, 'int8'));
   assertTrue(isa(A * x_Char, 'char'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function test_opNull_double(seed)
   n = randi(100);
   assertEqual( eye(n,'double'), double(opNull(n)) ) 
end

