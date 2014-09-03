function RxPower = Prop(TxPower, TxLoc, RxLocX, RxLocY)
    PL = 131.1 + 42.81*log10( sqrt((TxLoc(1) - RxLocX).^2 + (TxLoc(2) - RxLocY).^2)/1000 );
    RxPower = TxPower - PL;