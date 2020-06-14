using PyPlot
pygui(true)

KG = [1:500...]
KS = [1:500...]

PP = ones(length(KG), length(KS))
for i in range(1, length(KG), step = 1)
    for j in range(1, length(KS[i:end]), step = 1)
        PP[i, j] = KG[i]^KS[j]
    end
end

UL = ones(length(KG), length(KS))
for i in 1:length(KG)
    for j in 1:length(KS[i:end])
        UL[i, j] = factorial(big(KG[i]))/factorial(big(KG[i] - KS[i]))
    end
end

PP = [log(abs(pp)) for pp in PP]
UL = [log(abs(ul)) for ul in UL]

RR = zeros(length(KG), length(KS))
for i in 1:length(KG)
    for j in 1:length(KS)
        RR[i, j] = PP[i, j] / UL[i, j]
        if RR[i, j] == NaN
            RR[i, j] = 0
        end
    end
end

#[surf(KS[i], KG, RR) for i in 1:length(KS)]
surf(KS, KG, RR)
xlabel("K(S)")
ylabel("K(G)")
zlabel("Log Rel. Time")
