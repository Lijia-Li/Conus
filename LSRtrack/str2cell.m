%%%% NOTE: Courtesy of
%%%% http://mooring.ucsd.edu/software/matlab/mfiles/map/private/rmblank.m
%%%% DO NOT PUBLISH WITH THIS FILE UNLESS LISCENCING HAS BEEN CHECKED!

function v = str2cell(v)

% STR2CELL  Tries to Convert a String to CellArray
%
%  C = STR2CELL( String )
%
%  Try to build CellArray using EVAL(['{' String '}'])
%
%  In case of multiple Lines, separated by LF / CHAR(10)
%  a multiple Row CellArray for the nonempty Lines is returned.
%
% requires: RMBLANK
%

if ~chkstr(v,1) | isempty(v)
    return
end

ok = 1; try, v = eval(['{' v '}']); catch, ok = 0; end
if ok, return, end

%-----------------------------------------------------
% Remove Ending Blanks, NewLines etc.

    v = strrep(v,char( 0),' ');
    v = strrep(v,char( 9),' ');
    v = strrep(v,char(13),' ');

    v = rmblank(v,2);

ok = 1; try, v = eval(['{' v '}']); catch, ok = 0; end
if ok, return, end

%-----------------------------------------------------
% Remove ending "," ";"

   while any( v(end) == ',;' )
      if size(v,2) == 1
         v = '';
      else
         v = v( 1 : (end-1) );
         v = rmblank(v,2-i);
      end
      if isempty(v)
         break
      end
   end

ok = 1; try, v = eval(['{' v '}']); catch, ok = 0; end
if ok, return, end

%-----------------------------------------------------
% Remove surrounding "{}"

    if ~isempty(v)
        if v(1) == '{'
           if size(v,2) == 1
              v = '';
           else
              v = v( 2 : end );
              v = rmblank(v,2+i);
           end
        end
    end
    if ~isempty(v)
        if v(end) == '}'
           if size(v,2) == 1
              v = '';
           else
              v = v( 1 : (end-1) );
              v = rmblank(v,2-i);
           end
        end
    end

ok = 1; try, v = eval(['{' v '}']); catch, ok = 0; end
if ok, return, end

%-----------------------------------------------------
% Split by Lines

ok = ( v == 10 );

if ~any(ok)
    return
end

[i0,lg] = ind2grp(find(~ok));

ng = size(i0,1);

vv =  cell(ng,1);
nn = zeros(ng,1);

for ii = 1 : ng

    jj = ( 1 : lg(ii) ) + i0(ii) - 1;

    vv{ii} = rmblank(v(jj),2);

    if ~isempty(vv{ii})

        vv{ii} = str2cell(v(jj));

        nn(ii) = iscell(vv{ii});

        if nn(ii)
           if isempty(vv{ii})
              nn(ii) = -1;
           else
              vv{ii} =      vv{ii}(:)';
              nn(ii) = size(vv{ii},2);
           end
        end

    end

end

%-----------------------------------------------------
% Get Good Lines

ok = ~( nn == 0 );

if ~any(ok)
    return
end

n  = sum(ok);
ok = find(ok);

nn = nn(ok);
vv = vv(ok);

m  = max(nn);

if m == -1
   v = cell(n,0);
   return
end

v = cell(n,m);

for ii = 1 : n
    if ( nn(ii) == -1 )
        v(ii,1:nn(ii)) = vv{ii};
    end
end

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [i0,l] = ind2grp(ii);

% IND2GRP  Built StartIndex and Length from IndexVector
%
% [ StartIndex , GroupLength ] = IND2GRP( Index )
%

i0 = zeros(0,1);
l  = zeros(0,1);

if isempty(ii);
   return
end

ii = ii(:);
n  = size(ii,1);

i0 = cat( 1 , 1 , find( diff(ii,1,1) > 1 )+1 , n+1 );

l  = diff(i0,1,1);

i0 = ii(i0(1:end-1));

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );

   
   
function str = rmblank(str,dim,cc)

% RMBLANK  Remove Blanks, NewLines at Begin and End of CharacterArrays
%
% String = RMBLANK( CharArray )
%
% CharArray  2-dimensional CharacterArray
%
% further Options:
%
% String = RMBLANK( CharArray , DIM , CHAR )
%
%  
%  DIM  specifies Dimension to work, 
%       default: 2
%
%    A positive complex Part of DIM, to removes Blanks only from Start,
%    A negative complex Part of DIM, to removes Blanks only from End.
%       
%  CHAR specifies BlankCharacters to remove
%       default:  [ 160  32  13  10  9  0 ];  % [ NBSP Space CR LF TAB ZERO ]
%

  
msg = '';
 nl = char(10);

Nin = nargin;

if Nin < 1
  error('Not enough Input Arguments.')
else
  if ischar(str)
    str = double(str);
  end
  ok = isnumeric(str);
  if ok
    ok = all( ( mod(str(:),1) == 0 )  & ...
              ( str(:) >= 0 ) & isfinite(str(:))  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be a String or ASCII-Codes.'];
  end
  if size(str,1)*size(str,2) ~= prod(size(str))
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be 2-dimensional.'];
  end     
end

if Nin < 2
  dim = 2;
else
  if ~isnumeric(dim)
    msg = [ msg nl(1:(end*(~isempty(msg)))) ...
            'Input DIM must be numeric.' ];
  elseif ~isempty(dim)
    dim = dim(:);
    if ~all( ( real(dim) == 1 ) |  ( real(dim) == 2 ) )
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
             'Values for Input DIM must define 1. or 2. Dimension.' ];
    end
  end 
end

if Nin < 3
  cc = [ 160  32  13  10  9  0 ];  % [ NBSP  Space CR LF TAB ZERO ]
else
  if ischar(cc)
    cc = double(cc);
  end
  ok = isnumeric(cc);
  if ok & ~isempty(cc)
    cc = cc(:)';
    ok = all( ( mod(cc,1) == 0 )  & ...
              ( cc >= 0 ) & isfinite(cc)  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CHAR must be a String or ASCII-Codes.'];
  end
end

if ~isempty(msg)
  error(msg)
end


if isempty(str)
 str = '';
 return
end

if isempty(dim) | isempty(cc)
  str = double(str);
  return
end


  blank  = 0*str;

  for ii = cc
    blank = ( blank | ( str == ii ) );
  end

  si = size(str);

  for ii = 1 : size(dim,1)

    d = dim(ii);

    s = sign(imag(d));  % Remove from wich Side:  1  0  -1 
 
    d = real(d);

    jj = find( sum(blank,3-d) == si(3-d) );  % Columns with full Blanks

    if ~isempty(jj) 

         p  = [ 3-d  d ];
        str = permute(str,p);

         jj = jj(:)';
         nb = size(jj,2);

        %--------------------------------------------
        % Blank at Begin

        ind = ( 1 : nb );
        jj1 = find( ( ( jj == ind ) & ( s >= 0 ) ) );

        %--------------------------------------------
        % Blank at End

        ind = ind + si(d) - nb;
        jj2 = find( ( ( jj == ind ) & ( s <= 0 ) ) );

        %--------------------------------------------

        str(:,jj([jj1 jj2])) = [];

        str = permute(str,p);

    end
    
  end

  str = char(str);
