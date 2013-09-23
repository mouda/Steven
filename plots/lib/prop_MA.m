function [ RxPower ] = prop_MA( TxPower, TxLocX,TxLocY, RxLocX, RxLocY )
    PL = 131.1 + 42.68*log10( sqrt((TxLocX - RxLocX).^2 + (TxLocY - RxLocY).^2)/1000 );
    RxPower = TxPower - PL;
end

