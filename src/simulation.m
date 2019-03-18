% Define simulated SNR value
SNR_simulated = [ 5 , 10 , 15 , 20 , 25 ];
SNR_decoded = [];
BER_decoded = [];
result = [];

for i = 1 : length( SNR_simulated )
    signal_gen( SNR_simulated( i ) );
    [ SNR , BER ] = decode( SNR_simulated( i ) );
    result = [ result ; [ SNR , BER ] ];
end

for i = 1 : length( SNR_simulated )
    SNR_decoded = [ SNR_decoded , result( i , 1 ) ];
    BER_decoded = [ BER_decoded , result( i , 2 ) ];
end

% Step 5 - Plot SNR and BER
figure( 8 ); 
clf;

subplot( 2 , 1 , 1 );
bh = bar( SNR_simulated , SNR_decoded );
shading flat;
set( bh , 'FaceColor' , [ 0 , 0 , 1 ] );
%axis( [ min( SNR_simulated ) , max( SNR_simulated ) , 0 , max( SNR_decoded ) ] );
grid on;
title( 'SNR decoded' );
xlabel( 'Simulated SNR (dB)' );
ylabel( 'Actual SNR (dB)' );

subplot( 2 , 1 , 2 );
bh = bar( SNR_simulated , BER_decoded );
shading flat;
set( bh , 'FaceColor' , [ 0 , 0 , 1 ] );
%axis( [ min( SNR_simulated ) , max( SNR_simulated ) , 0 , max( BER_decoded ) + 5 ] );
grid on;
title( 'BER decoded' );
xlabel( 'Simulated SNR (dB)' );
ylabel( 'BER' );