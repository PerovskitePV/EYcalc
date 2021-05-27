% ------------- DOKUMENTATION OF THIS FUNCTION -------------
%
% #DESCRIPTION:           This function generates a sha1 hash to name the
%                         simulated optics files. Further a string "strline" 
%                         to write in a database and a generic data name 
%                         "DataNameOptics" is generated.
%
% #INPUT:                 input_string (string)
%
% #OUTPUT:                DataNameOptics (string)
%                         strline (string)
%                         hash (string)
%
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: -
%
% -----------------------------------------------------------

function hash = makehash( input_string )

    % Generate sha1 hash using .NET environement (only works with windows)
    sha1hasher = System.Security.Cryptography.SHA1Managed;
    sha1 = uint8(sha1hasher.ComputeHash(uint8(input_string)));
    hash = convertCharsToStrings(dec2hex(sha1));

end

