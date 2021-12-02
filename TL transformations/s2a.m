function abcd = s2a(s,Z_01,Z_02)
    
    S11 = squeeze(s(1,1,:));
    S12 = squeeze(s(1,2,:));
    S21 = squeeze(s(2,1,:));
    S22 = squeeze(s(2,2,:));
    
    den = 2.*S21 .* sqrt( real(Z_01) .* real(Z_02) );
    
    abcd(1,1,:) = (( conj(Z_01) + S11 .* Z_01 ) ./ ( 1 - S22 ) + S12 .* S21 .* Z_01 ) ./ den;
    abcd(1,2,:) = (( conj(Z_01) + S11 .* Z_01 ) .* ( conj(Z_02) + S22 .* Z_02 ) - S12 .* S21 .* Z_01 .* Z_02 ) ./ den;
    abcd(2,1,:) = (( 1 - S11 ) .* ( 1 - S22 ) - S12 .* S21 )./ den;
    abcd(2,2,:) = (( 1 - S11 ) .* ( conj(Z_02) + S22 .* Z_02 ) + S12 .* S21 .* Z_02 ) ./ den;
end