package gonudg

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/utils"
)

// GeometricFactors3D computes the metric elements for the local mappings of the elements
func (dg *DG3D) GeometricFactors3D() error {
	// Calculate geometric factors
	// xr = Dr*x, xs = Ds*x, xt = Dt*x
	xr := dg.Dr.Mul(dg.x)
	xs := dg.Ds.Mul(dg.x)
	xt := dg.Dt.Mul(dg.x)

	yr := dg.Dr.Mul(dg.y)
	ys := dg.Ds.Mul(dg.y)
	yt := dg.Dt.Mul(dg.y)

	zr := dg.Dr.Mul(dg.z)
	zs := dg.Ds.Mul(dg.z)
	zt := dg.Dt.Mul(dg.z)

	// Compute Jacobian determinant
	// J = xr*(ys*zt - zs*yt) - yr*(xs*zt - zs*xt) + zr*(xs*yt - ys*xt)
	dg.J = utils.NewMatrix(dg.Np, dg.K)
	for i := 0; i < dg.Np; i++ {
		for k := 0; k < dg.K; k++ {
			J := xr.At(i, k)*(ys.At(i, k)*zt.At(i, k)-zs.At(i, k)*yt.At(i, k)) -
				yr.At(i, k)*(xs.At(i, k)*zt.At(i, k)-zs.At(i, k)*xt.At(i, k)) +
				zr.At(i, k)*(xs.At(i, k)*yt.At(i, k)-ys.At(i, k)*xt.At(i, k))

			if J <= 0 {
				return fmt.Errorf("negative Jacobian at node %d, element %d: %f", i, k, J)
			}
			dg.J.Set(i, k, J)
		}
	}

	// Initialize metric term matrices
	dg.rx = utils.NewMatrix(dg.Np, dg.K)
	dg.ry = utils.NewMatrix(dg.Np, dg.K)
	dg.rz = utils.NewMatrix(dg.Np, dg.K)
	dg.sx = utils.NewMatrix(dg.Np, dg.K)
	dg.sy = utils.NewMatrix(dg.Np, dg.K)
	dg.sz = utils.NewMatrix(dg.Np, dg.K)
	dg.tx = utils.NewMatrix(dg.Np, dg.K)
	dg.ty = utils.NewMatrix(dg.Np, dg.K)
	dg.tz = utils.NewMatrix(dg.Np, dg.K)

	// Compute inverse metric terms
	for i := 0; i < dg.Np; i++ {
		for k := 0; k < dg.K; k++ {
			J := dg.J.At(i, k)

			// rx = (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz = (xs*yt - ys*xt)/J
			dg.rx.Set(i, k, (ys.At(i, k)*zt.At(i, k)-zs.At(i, k)*yt.At(i, k))/J)
			dg.ry.Set(i, k, -(xs.At(i, k)*zt.At(i, k)-zs.At(i, k)*xt.At(i, k))/J)
			dg.rz.Set(i, k, (xs.At(i, k)*yt.At(i, k)-ys.At(i, k)*xt.At(i, k))/J)

			// sx = -(yr*zt - zr*yt)/J, sy = (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J
			dg.sx.Set(i, k, -(yr.At(i, k)*zt.At(i, k)-zr.At(i, k)*yt.At(i, k))/J)
			dg.sy.Set(i, k, (xr.At(i, k)*zt.At(i, k)-zr.At(i, k)*xt.At(i, k))/J)
			dg.sz.Set(i, k, -(xr.At(i, k)*yt.At(i, k)-yr.At(i, k)*xt.At(i, k))/J)

			// tx = (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz = (xr*ys - yr*xs)/J
			dg.tx.Set(i, k, (yr.At(i, k)*zs.At(i, k)-zr.At(i, k)*ys.At(i, k))/J)
			dg.ty.Set(i, k, -(xr.At(i, k)*zs.At(i, k)-zr.At(i, k)*xs.At(i, k))/J)
			dg.tz.Set(i, k, (xr.At(i, k)*ys.At(i, k)-yr.At(i, k)*xs.At(i, k))/J)
		}
	}

	return nil
}

// Normals3D computes outward pointing normals at element faces and surface Jacobians
func (dg *DG3D) Normals3D() error {
	// First compute geometric factors
	if err := dg.GeometricFactors3D(); err != nil {
		return err
	}

	// Initialize normal and surface Jacobian matrices
	Nfp := dg.Nfp
	Nfaces := dg.Nfaces
	K := dg.K

	dg.nx = utils.NewMatrix(Nfp*Nfaces, K)
	dg.ny = utils.NewMatrix(Nfp*Nfaces, K)
	dg.nz = utils.NewMatrix(Nfp*Nfaces, K)
	dg.sJ = utils.NewMatrix(Nfp*Nfaces, K)

	// Build normals for each face
	// Face 1: t = -1, normal = -[tx, ty, tz]
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := dg.Fmask[0][i] // volume node index
			row := 0*Nfp + i      // face node index

			dg.nx.Set(row, k, -dg.tx.At(vid, k))
			dg.ny.Set(row, k, -dg.ty.At(vid, k))
			dg.nz.Set(row, k, -dg.tz.At(vid, k))
		}
	}

	// Face 2: s = -1, normal = -[sx, sy, sz]
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := dg.Fmask[1][i]
			row := 1*Nfp + i

			dg.nx.Set(row, k, -dg.sx.At(vid, k))
			dg.ny.Set(row, k, -dg.sy.At(vid, k))
			dg.nz.Set(row, k, -dg.sz.At(vid, k))
		}
	}

	// Face 3: r+s+t = -1, normal = [rx+sx+tx, ry+sy+ty, rz+sz+tz]
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := dg.Fmask[2][i]
			row := 2*Nfp + i

			dg.nx.Set(row, k, dg.rx.At(vid, k)+dg.sx.At(vid, k)+dg.tx.At(vid, k))
			dg.ny.Set(row, k, dg.ry.At(vid, k)+dg.sy.At(vid, k)+dg.ty.At(vid, k))
			dg.nz.Set(row, k, dg.rz.At(vid, k)+dg.sz.At(vid, k)+dg.tz.At(vid, k))
		}
	}

	// Face 4: r = -1, normal = -[rx, ry, rz]
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := dg.Fmask[3][i]
			row := 3*Nfp + i

			dg.nx.Set(row, k, -dg.rx.At(vid, k))
			dg.ny.Set(row, k, -dg.ry.At(vid, k))
			dg.nz.Set(row, k, -dg.rz.At(vid, k))
		}
	}

	// Normalize normals and compute surface Jacobian
	for i := 0; i < Nfp*Nfaces; i++ {
		for k := 0; k < K; k++ {
			nx := dg.nx.At(i, k)
			ny := dg.ny.At(i, k)
			nz := dg.nz.At(i, k)

			// Compute magnitude (surface Jacobian before scaling)
			mag := math.Sqrt(nx*nx + ny*ny + nz*nz)

			if mag < 1e-14 {
				return fmt.Errorf("zero normal magnitude at face node %d, element %d", i, k)
			}

			// Normalize to unit vector
			dg.nx.Set(i, k, nx/mag)
			dg.ny.Set(i, k, ny/mag)
			dg.nz.Set(i, k, nz/mag)

			// Store magnitude temporarily
			dg.sJ.Set(i, k, mag)
		}
	}

	// Scale surface Jacobian by volume Jacobian at face nodes
	// sJ = sJ * J(Fmask(:),:)
	for face := 0; face < Nfaces; face++ {
		for i := 0; i < Nfp; i++ {
			vid := dg.Fmask[face][i] // volume node index
			row := face*Nfp + i      // face node index

			for k := 0; k < K; k++ {
				dg.sJ.Set(row, k, dg.sJ.At(row, k)*dg.J.At(vid, k))
			}
		}
	}

	return nil
}
