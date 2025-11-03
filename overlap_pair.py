import numpy as np

def energy_theta(
    theta,
    h001, h100, h010, h011, h110, h020, h021, h120
):
    """
    Compute the energy corresponding to the long Fortran loop,
    but with u_p = cos(theta_p) and v_p = sin(theta_p).

    Parameters
    ----------
    theta : ndarray (na,)
        Angles for each site/orbital.
    h001, h100, h010 : ndarray (na,)
    h011, h110, h020 : ndarray (na, na)
    h021, h120       : ndarray (na, na, na)

    Returns
    -------
    E : float
        Total energy.
    """

    na = len(theta)
    u = np.cos(theta)
    v = np.sin(theta)

    # ---- tau0 ----
    tau0 = np.einsum('q,qp->p', v**2, h011)

    # ---- tau9 accumulates ----
    tau9 = 2 * tau0 * v

    # ---- tau1 ----
    tau1 = np.einsum('q,pq->p', v**2, h110)
    tau9 += 2 * tau1 * v

    # ---- tau2 ----
    tau2 = np.einsum('q,qp->p', (u**2)*(v**2), np.einsum('qqp->qp', h021))
    tau9 += 4 * tau2 * v

    # ---- tau3 ----
    tau3 = np.einsum('q,pqq->p', (u**2)*(v**2), h120)
    tau9 += 4 * tau3 * v

    # ---- tau4 ----
    tau4 = np.einsum('q,r, rqp -> p', v**2, v**2, h021)
    tau9 += 4 * tau4 * v

    # ---- tau5 ----
    tau5 = np.einsum('q,r, prq -> p', v**2, v**2, h120)
    tau9 += 4 * tau5 * v

    # ---- tau6 and tau7 ----
    tau6 = np.einsum('pq->pq', h021[:, :, :].diagonal(axis1=0, axis2=2)) + \
           np.einsum('ppq->pq', h120)
    tau7 = 4 * np.einsum('q,pq,q->p', np.ones(na), tau6, v**2)  # same as sum_q 4 v(q)^2 tau6(p,q)
    tau7 += np.diag(h011) + np.diag(h110)
    tau9 -= 2 * (v**3) * tau7

    # ---- tau8 ----
    diag021 = np.einsum('ppp->p', h021)
    diag120 = np.einsum('ppp->p', h120)
    tau8 = diag021 + diag120
    tau9 += 4 * (v**5) * tau8

    # ---- h001, h100 ----
    tau9 += h001 * v + h100 * v

    # ---- E partial ----
    E = np.dot(tau9, u)

    # ---- tau10 ----
    tau10 = h010 + 2*(u**2)*np.diag(h020)
    tau10 += 2 * np.einsum('q,pq->p', v**2, h020)
    E += np.dot(2*(v**2), tau10)

    # ---- tau11 ----
    tau11 = (v**3) * (diag021 + diag120)
    E -= np.dot(4*(u**3), tau11)

    return E
