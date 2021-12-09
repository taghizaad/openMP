function selectedBase = chooseBase(base1,base2,curTT)
base1Score = 0;
base2Score = 0;
    for i = 1:length(curTT)
        if(curTT(i) == base1(i))
            base1Score = base1Score + 1;
        end
        if (curTT(i) == base2(i))
            base2Score = base2Score + 1;
        end
    end
    if (base1Score >= base2Score) 
        selectedBase = base1;
    else 
        selectedBase = base2;
    end
end

