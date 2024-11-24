function f = fun_dispBardell(code,xi)

    code = upper(code);
    f = zeros(length(xi), 102);

    if code == "F"
        f(:,1)   =  0.500 - 0.750*xi + 0.250*xi.^3;
        f(:,2)   =  0.125 - 0.125*xi - 0.125*xi.^2 + 0.125*xi.^3;
        f(:,3)   =  0.500 + 0.750*xi - 0.250*xi.^3;
        f(:,4)   = -0.125 - 0.125*xi + 0.125*xi.^2 + 0.125*xi.^3;
        f(:,101) =  0.500 - 0.500*xi;
        f(:,102) =  0.500 + 0.500*xi;
        f(:,5)   =  0.125 - 0.250*xi.^2 + 0.125*xi.^4;
        f(:,6)   =  0.125*xi - 0.250*xi.^3 + 0.125*xi.^5;
        f(:,7)   =  0.14583333333333334*xi.^6 - 0.3125*xi.^4 + 0.1875*xi.^2 - 0.020833333333333332;
        f(:,8)   =  0.1875*xi.^7 - 0.4375*xi.^5 + 0.3125*xi.^3 - 0.0625*xi;
        f(:,9)   =  0.2578125*xi.^8 - 0.65625*xi.^6 + 0.546875*xi.^4 - 0.15625*xi.^2 + 0.0078125;
        f(:,10)  =  0.37239583333333331*xi.^9 - 1.03125*xi.^7 + 0.984375*xi.^5 - 0.36458333333333331*xi.^3 + 0.0390625*xi;
        f(:,11)  =  0.55859375*xi.^10 - 1.67578125*xi.^8 + 1.8046875*xi.^6 - 0.8203125*xi.^4 + 0.13671875*xi.^2 - 0.00390625;
        f(:,12)  =  0.86328125*xi.^11 - 2.79296875*xi.^9 + 3.3515625*xi.^7 - 1.8046875*xi.^5 + 0.41015625*xi.^3 - 0.02734375*xi;
        f(:,13)  =  1.3668619791666667*xi.^12 - 4.748046875*xi.^10 + 6.2841796875*xi.^8 - 3.91015625*xi.^6 + 1.1279296875*xi.^4 - 0.123046875*xi.^2 + 0.0022786458333333335;
        f(:,14)  =  2.2080078125*xi.^13 - 8.201171875*xi.^11 + 11.8701171875*xi.^9 - 8.37890625*xi.^7 + 2.9326171875*xi.^5 - 0.451171875*xi.^3 + 0.0205078125*xi;
        f(:,15)  =  3.62744140625*xi.^14 - 14.35205078125*xi.^12 + 22.55322265625*xi.^10 - 17.80517578125*xi.^8 + 7.33154296875*xi.^6 - 1.46630859375*xi.^4 + 0.11279296875*xi.^2 - 0.00146484375;
        f(:,16)  =  6.045735677083333*xi.^15 - 25.39208984375*xi.^13 + 43.05615234375*xi.^11 - 37.588704427083336*xi.^9 + 17.80517578125*xi.^7 - 4.39892578125*xi.^5 + 0.48876953125*xi.^3 - 0.01611328125*xi;
        f(:,17)  =  10.202178955078125*xi.^16 - 45.343017578125*xi.^14 + 82.5242919921875*xi.^12 - 78.936279296875*xi.^10 + 42.28729248046875*xi.^8 - 12.463623046875*xi.^6 + 1.8328857421875*xi.^4 - 0.104736328125*xi.^2 + 0.001007080078125;
        f(:,18)  =  17.403717041015625*xi.^17 - 81.617431640625*xi.^15 + 158.7005615234375*xi.^13 - 165.048583984375*xi.^11 + 98.67034912109375*xi.^9 - 33.829833984375*xi.^7 + 6.2318115234375*xi.^5 - 0.523681640625*xi.^3 + 0.013092041015625*xi;
        f(:,19)  =  29.973068237304688*xi.^18 - 147.93159484863281*xi.^16 + 306.06536865234375*xi.^14 - 343.85121663411456*xi.^12 + 226.94180297851563*xi.^10 - 88.803314208984375*xi.^8 + 19.73406982421875*xi.^6 - 2.22564697265625*xi.^4 + 0.0981903076171875*xi.^2 - 0.00072733561197916663;
        f(:,20)  =  52.058486938476563*xi.^19 - 269.75761413574219*xi.^17 + 591.72637939453125*xi.^15 - 714.15252685546875*xi.^13 + 515.77682495117188*xi.^11 - 226.94180297851563*xi.^9 + 59.20220947265625*xi.^7 - 8.45745849609375*xi.^5 + 0.5564117431640625*xi.^3 - 0.0109100341796875*xi;
        f(:,21)  =  91.102352142333984*xi.^20 - 494.55562591552734*xi.^18 + 1146.4698600769043*xi.^16 - 1479.3159484863281*xi.^14 + 1160.4978561401367*xi.^12 - 567.35450744628906*xi.^10 + 170.20635223388672*xi.^8 - 29.601104736328125*xi.^6 + 2.6429557800292969*xi.^4 - 0.09273529052734375*xi.^2 + 0.000545501708984375;
        f(:,22)  =  160.51366806030273*xi.^21 - 911.02352142333984*xi.^19 + 2225.500316619873*xi.^17 - 3057.2529602050781*xi.^15 + 2588.8029098510742*xi.^13 - 1392.5974273681641*xi.^11 + 472.79542287190753*xi.^9 - 97.260772705078125*xi.^7 + 11.100414276123047*xi.^5 - 0.58732350667317712*xi.^3 + 0.009273529052734375*xi;
        f(:,23)  =  284.54695701599121*xi.^22 - 1685.3935146331787*xi.^20 + 4327.3617267608643*xi.^18 - 6305.5842304229736*xi.^16 + 5732.3493003845215*xi.^14 - 3365.4437828063965*xi.^12 + 1276.5476417541504*xi.^10 - 303.93991470336914*xi.^8 + 42.55158805847168*xi.^6 - 3.0834484100341797*xi.^4 + 0.088098526000976563*xi.^2 - 0.0004215240478515625;
        f(:,24)  =  507.23587989807129*xi.^23 - 3130.0165271759033*xi.^21 + 8426.9675731658936*xi.^19 - 12982.085180282593*xi.^17 + 12611.168460845947*xi.^15 - 8025.2890205383301*xi.^13 + 3365.4437828063965*xi.^11 - 911.81974411010742*xi.^9 + 151.96995735168457*xi.^7 - 14.183862686157227*xi.^5 + 0.61668968200683594*xi.^3 - 0.0080089569091796875*xi;
        f(:,25)  =  908.79761815071106*xi.^24 - 5833.2126188278198*xi.^22 + 16432.586767673492*xi.^20 - 26685.39731502533*xi.^18 + 27586.93100810051*xi.^16 - 18916.752691268921*xi.^14 + 8694.0631055831909*xi.^12 - 2644.2772579193115*xi.^10 + 512.89860606193542*xi.^8 - 59.099427858988442*xi.^6 + 3.5459656715393066*xi.^4 - 0.084094047546386719*xi.^2 + 0.000333706537882487;
    elseif code == "DF"
        f(:,1)   = -0.750 + 0.750*xi.^2;
        f(:,2)   = -0.125 - 0.250*xi + 0.375*xi.^2;
        f(:,3)   = +0.750 - 0.750*xi.^2;
        f(:,4)   = -0.125 + 0.250*xi + 0.375*xi.^2;
        f(:,101) = -0.500;
        f(:,102) = +0.500;
        f(:,5)   = -0.500*xi + 0.500*xi.^3;
        f(:,6)   = 0.625*xi.^4 - 0.75*xi.^2 + 0.125;
        f(:,7)   = 0.875*xi.^5 - 1.25*xi.^3 + 0.375*xi;
        f(:,8)   = 1.3125*xi.^6 - 2.1875*xi.^4 + 0.9375*xi.^2 - 0.0625;
        f(:,9)   = 2.0625*xi.^7 - 3.9375*xi.^5 + 2.1875*xi.^3 - 0.3125*xi;
        f(:,10)  = 3.3515625*xi.^8 - 7.21875*xi.^6 + 4.921875*xi.^4 - 1.09375*xi.^2 + 0.0390625;
        f(:,11)  = 5.5859375*xi.^9 - 13.40625*xi.^7 + 10.828125*xi.^5 - 3.28125*xi.^3 + 0.2734375*xi;
        f(:,12)  = 9.49609375*xi.^10 - 25.13671875*xi.^8 + 23.4609375*xi.^6 - 9.0234375*xi.^4 + 1.23046875*xi.^2 - 0.02734375;
        f(:,13)  = 16.40234375*xi.^11 - 47.48046875*xi.^9 + 50.2734375*xi.^7 - 23.4609375*xi.^5 + 4.51171875*xi.^3 - 0.24609375*xi;
        f(:,14)  = 28.7041015625*xi.^12 - 90.212890625*xi.^10 + 106.8310546875*xi.^8 - 58.65234375*xi.^6 + 14.6630859375*xi.^4 - 1.353515625*xi.^2 + 0.0205078125;
        f(:,15)  = 50.7841796875*xi.^13 - 172.224609375*xi.^11 + 225.5322265625*xi.^9 - 142.44140625*xi.^7 + 43.9892578125*xi.^5 - 5.865234375*xi.^3 + 0.2255859375*xi;
        f(:,16)  = 90.68603515625*xi.^14 - 330.09716796875*xi.^12 + 473.61767578125*xi.^10 - 338.29833984375*xi.^8 + 124.63623046875*xi.^6 - 21.99462890625*xi.^4 + 1.46630859375*xi.^2 - 0.01611328125;
        f(:,17)  = 163.23486328125*xi.^15 - 634.80224609375*xi.^13 + 990.29150390625*xi.^11 - 789.36279296875*xi.^9 + 338.29833984375*xi.^7 - 74.78173828125*xi.^5 + 7.33154296875*xi.^3 - 0.20947265625*xi;
        f(:,18)  = 295.86318969726563*xi.^16 - 1224.261474609375*xi.^14 + 2063.1072998046875*xi.^12 - 1815.534423828125*xi.^10 + 888.03314208984375*xi.^8 - 236.808837890625*xi.^6 + 31.1590576171875*xi.^4 - 1.571044921875*xi.^2 + 0.013092041015625;
        f(:,19)  = 539.51522827148438*xi.^17 - 2366.905517578125*xi.^15 + 4284.9151611328125*xi.^13 - 4126.214599609375*xi.^11 + 2269.4180297851563*xi.^9 - 710.426513671875*xi.^7 + 118.4044189453125*xi.^5 - 8.902587890625*xi.^3 + 0.196380615234375*xi;
        f(:,20)  = 989.11125183105469*xi.^18 - 4585.8794403076172*xi.^16 + 8875.8956909179688*xi.^14 - 9283.9828491210938*xi.^12 + 5673.5450744628906*xi.^10 - 2042.4762268066406*xi.^8 + 414.41546630859375*xi.^6 - 42.28729248046875*xi.^4 + 1.6692352294921875*xi.^2 - 0.0109100341796875;
        f(:,21)  = 1822.0470428466797*xi.^19 - 8902.0012664794922*xi.^17 + 18343.517761230469*xi.^15 - 20710.423278808594*xi.^13 + 13925.974273681641*xi.^11 - 5673.5450744628906*xi.^9 + 1361.6508178710938*xi.^7 - 177.60662841796875*xi.^5 + 10.571823120117188*xi.^3 - 0.1854705810546875*xi;
        f(:,22)  = 3370.7870292663574*xi.^20 - 17309.446907043457*xi.^18 + 37833.505382537842*xi.^16 - 45858.794403076172*xi.^14 + 33654.437828063965*xi.^12 - 15318.571701049805*xi.^10 + 4255.158805847168*xi.^8 - 680.82540893554688*xi.^6 + 55.502071380615234*xi.^4 - 1.7619705200195313*xi.^2 + 0.009273529052734375;
        f(:,23)  = 6260.0330543518066*xi.^21 - 33707.870292663574*xi.^19 + 77892.511081695557*xi.^17 - 100889.34768676758*xi.^15 + 80252.890205383301*xi.^13 - 40385.325393676758*xi.^11 + 12765.476417541504*xi.^9 - 2431.5193176269531*xi.^7 + 255.30952835083008*xi.^5 - 12.333793640136719*xi.^3 + 0.17619705200195313*xi;
        f(:,24)  = 11666.42523765564*xi.^22 - 65730.34707069397*xi.^20 + 160112.38389015198*xi.^18 - 220695.44806480408*xi.^16 + 189167.52691268921*xi.^14 - 104328.75726699829*xi.^12 + 37019.881610870361*xi.^10 - 8206.3776969909668*xi.^8 + 1063.789701461792*xi.^6 - 70.919313430786133*xi.^4 + 1.8500690460205078*xi.^2 - 0.0080089569091796875;
        f(:,25)  = 21811.142835617065*xi.^23 - 128330.67761421204*xi.^21 + 328651.73535346985*xi.^19 - 480337.15167045593*xi.^17 + 441390.89612960815*xi.^15 - 264834.53767776489*xi.^13 + 104328.75726699829*xi.^11 - 26442.772579193115*xi.^9 + 4103.1888484954834*xi.^7 - 354.59656715393066*xi.^5 + 14.183862686157227*xi.^3 - 0.16818809509277344*xi;
    elseif code == "D2F"
        f(:,1)   = +1.500*xi;
        f(:,2)   = -0.250 + 0.750*xi;
        f(:,3)   = -1.500*xi;
        f(:,4)   = +0.250 + 0.750*xi;
        f(:,101) = 0;
        f(:,102) = 0;
        f(:,5)   = 1.5*xi.^2- 0.5;
        f(:,6)   = xi.*(2.5*xi.^2- 1.5);
        f(:,7)   = 4.375*xi.^4- 3.75*xi.^2+ 0.375;
        f(:,8)   = xi.*(7.875*xi.^4- 8.75*xi.^2+ 1.875);
        f(:,9)   = 14.4375*xi.^6- 19.6875*xi.^4+ 6.5625*xi.^2- 0.3125;
        f(:,10)  = xi.*(26.8125*xi.^6- 43.3125*xi.^4+ 19.6875*xi.^2- 2.1875);
        f(:,11)  = 50.2734375*xi.^8- 93.84375*xi.^6+ 54.140625*xi.^4- 9.84375*xi.^2+ 0.2734375;
        f(:,12)  = xi.*(94.9609375*xi.^8- 201.09375*xi.^6+ 140.765625*xi.^4- 36.09375*xi.^2+ 2.4609375);
        f(:,13)  = 180.42578125*xi.^10- 427.32421875*xi.^8+ 351.9140625*xi.^6- 117.3046875*xi.^4+ 13.53515625*xi.^2- 0.24609375;
        f(:,14)  = xi.*(344.44921875*xi.^10- 902.12890625*xi.^8+ 854.6484375*xi.^6- 351.9140625*xi.^4+ 58.65234375*xi.^2- 2.70703125);
        f(:,15)  = 660.1943359375*xi.^12- 1894.470703125*xi.^10+ 2029.7900390625*xi.^8- 997.08984375*xi.^6+ 219.9462890625*xi.^4- 17.595703125*xi.^2+ 0.2255859375;
        f(:,16)  = xi.*(1269.6044921875*xi.^12- 3961.166015625*xi.^10+ 4736.1767578125*xi.^8- 2706.38671875*xi.^6+ 747.8173828125*xi.^4- 87.978515625*xi.^2+ 2.9326171875);
        f(:,17)  = 2448.52294921875*xi.^14- 8252.42919921875*xi.^12+ 10893.20654296875*xi.^10- 7104.26513671875*xi.^8+ 2368.08837890625*xi.^6- 373.90869140625*xi.^4+ 21.99462890625*xi.^2- 0.20947265625;
        f(:,18)  = xi.*(4733.81103515625*xi.^14- 17139.66064453125*xi.^12+ 24757.28759765625*xi.^10- 18155.34423828125*xi.^8+ 7104.26513671875*xi.^6- 1420.85302734375*xi.^4+ 124.63623046875*xi.^2- 3.14208984375);
        f(:,19)  = 9171.7588806152344*xi.^16- 35503.582763671875*xi.^14+ 55703.897094726563*xi.^12- 45388.360595703125*xi.^10+ 20424.762268066406*xi.^8- 4972.985595703125*xi.^6+ 592.0220947265625*xi.^4- 26.707763671875*xi.^2+ 0.196380615234375;
        f(:,20)  = xi.*(17804.002532958984*xi.^16- 73374.071044921875*xi.^14+ 124262.53967285156*xi.^12- 111407.79418945313*xi.^10+ 56735.450744628906*xi.^8- 16339.809814453125*xi.^6+ 2486.4927978515625*xi.^4- 169.149169921875*xi.^2+ 3.338470458984375);
        f(:,21)  = 34618.893814086914*xi.^18- 151334.02153015137*xi.^16+ 275152.76641845703*xi.^14- 269235.50262451172*xi.^12+ 153185.71701049805*xi.^10- 51061.905670166016*xi.^8+ 9531.5557250976563*xi.^6- 888.03314208984375*xi.^4+ 31.715469360351563*xi.^2- 0.1854705810546875;
        f(:,22)  = xi.*(67415.740585327148*xi.^18- 311570.04432678223*xi.^16+ 605336.08612060547*xi.^14- 642023.12164306641*xi.^12+ 403853.25393676758*xi.^10- 153185.71701049805*xi.^8+ 34041.270446777344*xi.^6- 4084.9524536132813*xi.^4+ 222.00828552246094*xi.^2- 3.5239410400390625);
        f(:,23)  = 131460.69414138794*xi.^20- 640449.53556060791*xi.^18+ 1324172.6883888245*xi.^16- 1513340.2153015137*xi.^14+ 1043287.5726699829*xi.^12- 444238.57933044434*xi.^10+ 114889.28775787354*xi.^8- 17020.635223388672*xi.^6+ 1276.5476417541504*xi.^4- 37.001380920410156*xi.^2+ 0.17619705200195313;
        f(:,24)  = xi.*(256661.35522842407*xi.^20- 1314606.9414138794*xi.^18+ 2882022.9100227356*xi.^16- 3531127.1690368652*xi.^14+ 2648345.3767776489*xi.^12- 1251945.0872039795*xi.^10+ 370198.81610870361*xi.^8- 65651.021575927734*xi.^6+ 6382.738208770752*xi.^4- 283.67725372314453*xi.^2+ 3.7001380920410156);
        f(:,25)  = 501656.2852191925*xi.^22- 2694944.2298984528*xi.^20+ 6244382.9717159271*xi.^18- 8165731.5783977509*xi.^16+ 6620863.4419441223*xi.^14- 3442848.9898109436*xi.^12+ 1147616.3299369812*xi.^10- 237984.95321273804*xi.^8+ 28722.321939468384*xi.^6- 1772.9828357696533*xi.^4+ 42.55158805847168*xi.^2- 0.16818809509277344;
    end

end