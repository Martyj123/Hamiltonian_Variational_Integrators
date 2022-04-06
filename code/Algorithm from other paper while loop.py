while t < 20:
    #Paper based on Pihajoki algorithm, not mine - missing some equations
    q1half = q11 + dhdp1(q11t, q21t, p11, p21) * dt
    p1halft = (p11t + dhdq1(q11t, q21t, p11, p21) * dt)
    q2half = q21 + dhdp2(q11t, q21t, p11, p21) * dt
    p2halft = (p21t + dhdq2(q11t, q21t, p11, p21) * dt)

    q1halft = q11t + dhdp1(q1half, q2half, p1halft, p2halft) * dt
    p1half = (p11 + dhdq1(q1half, q2half, p1halft, p2halft) * dt)
    p2half = (p21 + dhdq2(q1half, q2half, p1halft, p2halft) * dt)
    q2halft = q21t + dhdp2(q1half, q2half, p1halft, p2halft) * dt

    p1n1end = p1half + dhdq1(q1half, q2half, p1halft, p2halft) * dt
    q1n1endt = q1halft + dhdp1(q1half, q2half, p1halft, p2halft) * dt
    p2n1end = p2half + dhdq2(q1half, q2half, p1halft, p2halft) * dt
    q2n1endt = q2halft + dhdp2(q1half, q2half, p1halft, p2halft) * dt
    
    p1n1endt = p1halft + dhdq1(q1n1endt, q2n1endt, p1n1end, p2n1end) * dt
    q1n1end = q1half + dhdp1(q1n1endt, q2n1endt, p1n1end, p2n1end) * dt
    p2n1endt = p2halft + dhdq2(q1n1endt, q2n1endt, p1n1end, p2n1end) * dt
    q2n1end = q2half + dhdp2(q1n1endt, q2n1endt, p1n1end, p2n1end) * dt
    
    #Calculate/store positions for plotting
    posx1= l1 * np.sin(q1n1end)
    posy1= -l1 * np.cos(q1n1end)
    posx2= posx1 + l2 * np.sin(q2n1end)
    posy2= posy1 - l2 * np.cos(q2n1end)
    posxvals1.append(posx1)
    posyvals1.append(posy1)
    posxvals2.append(posx2)
    posyvals2.append(posy2)
    #Update values for next time step
    q11 = q1n1end
    q21 = q2n1end
    p11 = p1n1end
    p21 = p2n1end
    q11t = q1n1end
    q21t = q2n1end
    p11t = p1n1end
    p21t = p2n1end
    n+=1
    t+=dt;
    tvals.append(t)
    E=hamilt(q1n1end,q2n1end,p1n1end,p2n1end)
    evals.append(E)