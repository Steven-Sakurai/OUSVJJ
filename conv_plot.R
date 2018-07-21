lb = c(-0.5, 0.001, 1e-6, 0.001, -0.95, 0, 0, -0.95, -1, 0, 0, 0)
ub = c(1, 100, 10, 1, 0.95, 10, 0, 0, 0, 1, 1, 1)

converge.plot.n <- function(est.par, i) {
    n = 12
    name.list = c("mu", "kappa", "theta", "alpha", "rho", "l0", "l1", "rho_j", "mu_s", "sigma_s", "mu_v", "sigma_v")    
    f.test = c()
    par.test = c()
    
    add = rep(0, n)
    if(abs(est.par[i]) >= 1e-5) {
        add[i] = est.par[i]
        add = sapply(add, function(x) x*seq(-0.08, 0.08, 0.01))
    } 
    else if(lb[i] == ub[i]) {
        add = sapply(add, function(x) x*seq(-0.08, 0.08, 0.01))
    }
    else {
        add[i] = 0.001
        add = sapply(add, function(x) x*seq(-8, 8, 1))
    }
    
    for(j in 1:(dim(add)[1])) {
        par.test = c(par.test, est.par[i] + add[j, i])
        f.test = c(f.test, f(est.par + add[j, ]))
    }
    df1 = data.frame(par.test, f.test, size = c(rep(0, 8), 2, rep(0, 8)))
    p1 = ggplot(df1, aes(x = par.test, y = f.test)) + geom_line() + geom_point(aes(color = "red", size = size))
    p1 = p1 + xlab(name.list[i]) + ylab("nll") + scale_size(guide = "none") + scale_color_discrete(guide = "none")
    return(p1)
}

converge.plot <- function(est.par) {
    p1 = converge.plot.n(est.par, 1)
    p2 = converge.plot.n(est.par, 2)
    p3 = converge.plot.n(est.par, 3)
    p4 = converge.plot.n(est.par, 4)
    p5 = converge.plot.n(est.par, 5)
    p6 = converge.plot.n(est.par, 6)
    p7 = converge.plot.n(est.par, 7)
    p8 = converge.plot.n(est.par, 8)
    p9 = converge.plot.n(est.par, 9)
    p10 = converge.plot.n(est.par, 10)
    p11 = converge.plot.n(est.par, 11)
    p12 = converge.plot.n(est.par, 12)
    return(grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, ncol = 3))
}

check.converge.n <- function(f, par, n) {
    f0 = f(par)
    if(abs(par[n]) <= 1e-5) {
        par1 = par; par1[n] = 0.001; f1 = f(par1)
        par3 = par; par3[n] = 0.002; f3 = f(par3)

        par2 = par; par2[n] = -0.001; f2 = f(par2)
        par4 = par; par4[n] = -0.002; f4 = f(par4)
    }
    else {
        par1 = par; par1[n] = par1[n] + abs(par1[n]*2e-2); f1 = f(par1)
        par3 = par; par3[n] = par3[n] + abs(par3[n]*4e-2); f3 = f(par3)

        par2 = par; par2[n] = par2[n] - abs(par2[n]*2e-2); f2 = f(par2)
        par4 = par; par4[n] = par4[n] - abs(par4[n]*4e-2); f4 = f(par4)
    }
    f
    cond_lst = c(f1 >= f0, f3 >= f0, f2 >= f0, f4 >= f0)
    cond_lst = cond_lst[ !is.nan(cond_lst) ]
    return(all(cond_lst))
}
check.converge <- function(f, sol) {
    a1 = check.converge.n(f, sol, 1)
    a2 = check.converge.n(f, sol, 2)
    a3 = check.converge.n(f, sol, 3)
    a4 = check.converge.n(f, sol, 4)
    a5 = check.converge.n(f, sol, 5)
    a6 = check.converge.n(f, sol, 6)
    a7 = check.converge.n(f, sol, 7)
    a8 = check.converge.n(f, sol, 8)
    a9 = check.converge.n(f, sol, 9)
    a10 = check.converge.n(f, sol, 10)
    a11 = check.converge.n(f, sol, 11)
    a12 = check.converge.n(f, sol, 12)
    return(c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12))
}

check.mono.n <- function(f, par, n) {
    lower = lb[n]
    upper = ub[n]
    c1 = (par[n] < lower + 0.001 + abs(lower*0.01))
    c2 = (par[n] > upper - 0.001 - abs(upper*0.01))

    f0 = f(par)
    if(abs(par[n]) <= 1e-5) {
    	par1 = par; par1[n] = 0.001; f1 = f(par1)
	    par3 = par; par3[n] = 0.002; f3 = f(par3)

	    par2 = par; par2[n] = -0.001; f2 = f(par2)
	    par4 = par; par4[n] = -0.002; f4 = f(par4)
    }
    else {
    	par1 = par; par1[n] = par1[n] + abs(par1[n]*1e-2); f1 = f(par1)
	    par3 = par; par3[n] = par3[n] + abs(par3[n]*2e-2); f3 = f(par3)

	    par2 = par; par2[n] = par2[n] - abs(par2[n]*1e-2); f2 = f(par2)
	    par4 = par; par4[n] = par4[n] - abs(par4[n]*2e-2); f4 = f(par4)
    }
    cond1 = TRUE
    cond2 = TRUE
    if(!( is.nan(f1) || is.nan(f3)) ) {
        cond1 = cond1 && f3 >= f1
        cond2 = cond2 && f3 <= f1
    }
    if(!( is.nan(f2) || is.nan(f4)) ) {
        cond1 = cond1 && f2 >= f4
        cond2 = cond2 && f2 <= f4
    }
    if(!( is.nan(f1)) ) {
        cond1 = cond1 && f1 >= f0
        cond2 = cond2 && f1 <= f0
    }
    if(!( is.nan(f2)) ) {
        cond1 = cond1 && f2 <= f0
        cond2 = cond2 && f2 >= f0
    }
    return( (cond1 && c1) || (cond2 && c2) )
}

check.mono <- function(f, sol) {
    a1 = check.mono.n(f, sol, 1)
    a2 = check.mono.n(f, sol, 2)
    a3 = check.mono.n(f, sol, 3)
    a4 = check.mono.n(f, sol, 4)
    a5 = check.mono.n(f, sol, 5)
    a6 = check.mono.n(f, sol, 6)
    a7 = check.mono.n(f, sol, 7)
    a8 = check.mono.n(f, sol, 8)
    a9 = check.mono.n(f, sol, 9)
    a10 = check.mono.n(f, sol, 10)
    a11 = check.mono.n(f, sol, 11)
    a12 = check.mono.n(f, sol, 12)
    ret = c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12)
    return(ret)
}

