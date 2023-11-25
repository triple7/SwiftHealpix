//
//  Healpix.swift
//  swift-healpix
//
//  Created by Yuma decaux on 2/3/2022.
//

import Foundation
import simd
import BigInt

public typealias HP = Healpix

public final class Healpix {
    private let NSIDELIST:[Int] = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216]
    private let NS_MAX = 16777216
static let ORDER_MAX = 24
    static let PI2 = 2 * Double.pi
    static let PI = Double.pi
static let PI_2 = Double.pi / 2
    static let PI_4 = Double.pi / 4
    static let PI_8 = Double.pi / 8;

    private let nside:Int!
    
    public init(_ nside: Int) {
        self.nside = nside
    }
    
    public final class func order2Nside(_ order: Int)->Int {
        return 1 << order
    }

    public final class func nside2Order(_ nside: Int)->Int {
        return Int(log2(Float(nside)))
    }
    
    public final class func nside2Npix(_ nside: Int)->Int {
        return 12*nside*nside
    }
    
    public final class func vec2Pix_nest(_ nside: Int, _ v: simd_double3)->Int {
        let za = Healpix.vec2Za(v)
        return Healpix.za2Pix_nest(nside, za)
    }
    
    public final class func vec2Pix_ring(_ nside: Int, _ v: simd_double3)->Int {
        let za = Healpix.vec2Za(v)
        return Healpix.nest2Ring(nside, Healpix.za2Pix_nest(nside, za))
    }
    
    public final class func ang2Pix_nest(_ nside: Int, _ theta: Double, _ phi: Double)->Int {
        let z = cos(theta)
        return Healpix.za2Pix_nest(nside, simd_double2([z, phi]))
    }
    
    public final class func ang2Pix_ring(_ nside: Int, _ theta: Double, _ phi: Double)->Int {
        let z = cos(theta)
        return Healpix.nest2Ring(nside, Healpix.za2Pix_nest(nside, simd_double2([z, phi])))
    }
    
    public final class func nest2Ring(_ nside: Int, _ ipix: Int)->Int {
        let fxy = Healpix.nest2Fxy(nside, ipix)
        return Healpix.fxy2Ring(nside, fxy)
    }
    
    public final class func ring2Nest(_ nside: Int, _ ipix: Int)->Int {
        if nside == 1 {
            return ipix
        }
        let fxy = Healpix.ring2Fxy(nside, ipix)
        return Healpix.fxy2Nest(nside, fxy)
    }

    public final class func ring2Fxy(_ nside: Int, _ ipix: Int)->simd_int3 {
        let polarLimit = 2*nside*(nside - 1)
        if ipix < polarLimit {
            let ipix = Double(ipix)
            let nside = Double(nside)
            let i = floor((sqrt(1 + 2*ipix) + 1) / 2)
            let j = ipix - Double(2*i*(i - 1))
            let f = floor(j/i)
            let k = j.remainder(dividingBy: i)
            let x = nside - i + k
            let y = nside - 1 - k
            return simd_int3(Int32(f), Int32(x), Int32(y))
        }
        if ipix < polarLimit + 8*nside*nside {
            let nside = Double(nside)
            let k = Double(ipix - polarLimit)
            let ring = 4*nside
            let i = nside - floor(k / ring)
            let s:Double = (i.remainder(dividingBy: 2) == 0) ? 1 : 0
            let j = 2*k.remainder(dividingBy: ring) + s
            let jj = j - 4*nside
            let ii = i + 5*nside - 1
            let pp = (ii + jj) / 2
            let qq = (ii - jj) / 2
            let PP = floor(pp / nside)
            let QQ = floor(qq / nside)
            let V = 5 - (PP + QQ)
            let H = PP - QQ + 4
            let f = 4*Int(V) + (Int(H) >> 1) % 4
            let x = pp.remainder(dividingBy: nside)
            let y = qq.remainder(dividingBy: nside)
            return simd_int3(Int32(f), Int32(x), Int32(y))
        } else {
            let nside = Double(nside)
            let ipix = Double(ipix)
            let p = 12*nside*nside - ipix - 1
            let i = floor((sqrt(1 + 2 * p) + 1) / 2)
            let j = p - 2 * i * (i - 1)
            let f = 11 - floor(j / i)
            let k = j.remainder(dividingBy: i)
            let x = i - k - 1;
            let y = k
            return simd_int3(Int32(f), Int32(x), Int32(y))
        }
    }

    public final class func pix2Vec_nest(_ nside: Int, _ ipix: Int)->simd_double3 {
        let fxy = Healpix.nest2Fxy(nside, ipix)
        let tu = Healpix.fxy2Tu(nside, fxy)
        let za = Healpix.tu2Za(tu)
        return Healpix.za2Vec(za)
    }
    
    public final class func pix2Ang_nest(_ nside: Int, _ ipix: Int)->simd_double2 {
        let fxy = Healpix.nest2Fxy(nside, ipix)
        let tu = Healpix.fxy2Tu(nside, fxy)
        let za = Healpix.tu2Za(tu)
        return simd_double2([acos(za.x), za.y])
    }
    
    public final class func pix2Vec_ring(_ nside: Int, _ ipix: Int)->simd_double3 {
        return Healpix.pix2Vec_nest(nside, Healpix.ring2Nest(nside, ipix))
    }
    
    public final class func pix2Ang_ring(_ nside: Int, _ ipix: Int)->simd_double2 {
        return Healpix.pix2Ang_nest(nside, Healpix.ring2Nest(nside, ipix))
    }

    final func queryDisc(_ v: simd_double3, _  radius: Double, cb: (Int)->[Int])->[Int] {
        let nside = self.nside!
        if radius > Healpix.PI_2 {
            fatalError("Radius must be under PI/2")
        }
        let rad = Healpix.max_pixrad(nside)
        let d = Healpix.PI_4 / Double(nside)
        let za = Healpix.vec2Za(v)
        let sin_t = sqrt(1 - za.x * za.x)
        let cos_r = cos(radius)
        let sin_r = sin(radius)
        let z1 = za.x*cos_r + sin_t*sin_r // cos(theta - r)
        let z2 = za.x*cos_r - sin_t*sin_r
        let u1 = Healpix.za2Tu(simd_double2([z1, 0])).y
        let u2 = Healpix.za2Tu(simd_double2([z2, 0])).y
        let cover_north_pole = sin_t*cos_r - za.x*sin_r < 0
        let cover_south_pole = sin_t * cos_r + za.x * sin_r < 0
        var i1 = floor((Healpix.PI_2 - u1) / d)
        var i2 = floor((Healpix.PI_2 - u2) / d + 1)
        if cover_north_pole {
            i1 += 1
            for i in 1 ... Int(i1) {
                Healpix.walk_ring(nside, i, cb)
        }
            i1 += 1
        }
        if i1 == 0 {
            Healpix.walk_ring(nside, 1, cb)
            i1 = 2
        }
        if cover_south_pole {
            i2 -= 1
            for i in Int(i2) ... 4*nside-1 {
                Healpix.walk_ring(nside, i, cb)
        }
            i2 -= 1
        }
        if Int(i2) == 4*nside {
            Healpix.walk_ring(nside, 4 * nside - 1, cb)
            i2 = 4*Double(nside) - 2
        }
        let theta = acos(za.x)
        var output = [Int]()
        for i in Int(i1) ... Int(i2) {
            Healpix.walk_ring_around(nside, i, za.y, theta, radius + rad, { ipix in
                if (Healpix.angle(Healpix.pix2Vec_nest(nside, ipix), v) <= radius + rad) {
                    output += cb(ipix)
                }
                return output
            })
        }
        return output
    }

    public final class func max_pixrad(_ nside: Int)->Double {
        let nside = Double(nside)
        let unit = Healpix.PI_4 / nside
        return Healpix.angle(
            Healpix.tu2Vec(simd_double2(unit, nside*unit)),
            Healpix.tu2Vec(simd_double2(unit, (nside + 1)*unit))
        )
    }
    
    public final class func angle(_ a: simd_double3, _ b: simd_double3)->Double {
        return 2*asin(sqrt(simd_distance(a, b)) / 2)
    }
    
    public final class func tu2Vec(_ tu: simd_double2)->simd_double3 {
        let za = Healpix.tu2Za(tu)
        return Healpix.za2Vec(za)
    }
    
    public final class func walk_ring_around(_ nside: Int, _ i: Int, _ a0: Double, _ theta: Double, _ r: Double, _ cb: (Int)->[Int]) {
        if theta < r || theta + r > Double.pi {
            return walk_ring(nside, i, cb)
        }
        let nside = Double(nside)
        let u = Healpix.PI_4*(2 - Double(i) / nside)
        let z = Healpix.tu2Za(simd_double2(PI_4, u)).x
        let st = sin(theta)
        let ct = cos(theta)
        let sr = sin(r)
        let cr = cos(r)
        let asqr = z - ct * cr
        let bsqr = st*st * sr*sr
        let w = atan2(
            sqrt(-asqr / bsqr + 1) * sr,
            (-z * ct + cr) / st)
        let nsideInt = Int(nside)

        if w >= Double.pi {
            return Healpix.walk_ring(nsideInt, i, cb)
        }
        let t1 = Healpix.center_t(nsideInt, i, Healpix.za2Tu(simd_double2(z, Double(Healpix.wrap(a0 - w,Healpix.PI2)))).x)
        let t2 = Healpix.center_t(nsideInt, i, Healpix.za2Tu(simd_double2(z, Double(Healpix.wrap(a0 + w,Healpix.PI2)))).x)
        let begin = Healpix.tu2Fxy(nsideInt, simd_double2(t1, u))
        let end = Healpix.write_next_pixel(nsideInt, Healpix.tu2Fxy(nsideInt, simd_double2(t2, u)))
        var s = begin
        while !simd_equal(s, end) {
            var _ = cb(Healpix.fxy2Nest(nsideInt, s))
            s = Healpix.write_next_pixel(nsideInt, s)
        }
    }

    public final class func walk_ring(_ nside: Int, _ i: Int, _ cb: (Int)->[Int]) {
            let u = Healpix.PI_4 * (2 - Double(i) / Double(nside))
            let t = Healpix.PI_4 * (1 + Double(1 - i % 2) / Double(nside))
            let begin = Healpix.tu2Fxy(nside, simd_double2(t, u))
            var s = begin
            repeat {
                var _ = cb(Healpix.fxy2Nest(nside, s))
                s = Healpix.write_next_pixel(nside, s)
            } while !simd_equal(s, begin)
        }

    public final class func center_t(_ nside: Int, _ i: Int, _ t: Double)->Double {
        var t = Double(t)
            let d = Healpix.PI_4 / Double(nside)
            t /= d
            let t1 = (((Int(t) + i % 2) >> 1) << 1) + 1 - i % 2
            return Double(t1) * d
        }

    
    public final class func write_next_pixel(_ nside: Int, _ fxy: simd_int3)->simd_int3 {
        var fxy = fxy
        fxy.y += 1
        if fxy.y == nside {
            switch fxy.y/4 {
            case 0:
                fxy.x = (fxy.x + 1) % 4
                fxy.y = fxy.z
                fxy.z = Int32(nside)
            case 1:
                fxy.x = fxy.x - 4
                fxy.y = 0
            case 2:
                fxy.x = 4 + (fxy.x + 1) % 4
                fxy.y = 0
            default:
                break
            }
        }
        fxy.z -= 1
        if fxy.z == -1 {
            switch fxy.x / 4 {
            case 0:
                fxy.x = 4 + (fxy.x + 1) % 4
                fxy.z = Int32(nside) - 1
            case 1:
                fxy.z = Int32(nside) - 1
            case 2:
                fxy.x = 8 + (fxy.x + 1) % 4
                fxy.z = fxy.y - 1
                fxy.y = 0
            default:
                break
            }
        }
        return fxy
    }
    
    public final class func corners_nest(_ ipix: Int, _ nside: Int)->[simd_int3] {
        let fxy = Healpix.nest2Fxy(nside, ipix)
        let tu = Healpix.fxy2Tu(nside, fxy)
        let d = Healpix.PI_4 / Double(nside)
        var xyzs = [simd_double3]()
        for (tt, uu) in zip([0, d, 0, d], [d, 0, d, 0]) {
            let za = Healpix.tu2Za( simd_double2(tu.x + tt, tu.y + uu))
            xyzs.append(Healpix.za2Vec(za))
        }
        return xyzs.map{simd_int3(Int32($0.x), Int32($0.y), Int32($0.z))}
    }
             

    public final class func corners_ring(_ nside: Int, _ ipix: Int)->[simd_int3] {
            return Healpix.corners_nest(nside, Healpix.ring2Nest(nside, ipix))
        }
             
             public final class func nside2PixArea(_ nside: Int)->Double {
            return Double.pi / Double(3*nside*nside)
        }
             
    public final class func pixArea2Nside(_ pixArea: Double)->Int {
   return Int(sqrt(Double.pi / Double(3*pixArea)))
}
             public final class func nside2Resol(_ nside: Int)->Double {
            return sqrt(Double.pi/3) / Double(nside)
        }
             
             public final class func pixCoord2Vec_nest(_ nside: Int, _ ipix: Int, _ ne: Int, _ nw: Int)->simd_int3 {
            let fxy = Healpix.nest2Fxy(nside, ipix)
            let tu = Healpix.fxy2Tu(nside, fxy)
            let d = Healpix.PI_4 / Double(nside)
            let za = Healpix.tu2Za(simd_double2(tu.x + d*Double(ne - nw), tu.y + d*Double(ne + nw - 1)))
            let v = Healpix.za2Vec(za)
                 return simd_int3(Int32(v.x), Int32(v.y), Int32(v.z))
        }

             public final class func pixCoord2Vec_ring(_ nside: Int, _ ipix: Int, _ ne: Int, _ nw: Int)->simd_int3 {
            return Healpix.pixCoord2Vec_nest(nside, Healpix.ring2Nest(nside, ipix), ne, nw)
        }
             
             public final class func za2Pix_nest(_ nside: Int, _ za: simd_double2)->Int {
            let tu = Healpix.za2Tu(za)
            let fxy = Healpix.tu2Fxy(nside, tu)
            return Healpix.fxy2Nest(nside, fxy)
        }
             
             
             public final class func tu2Fxy(_ nside: Int, _ tu: simd_double2)->simd_int3 {
            let fpq = Healpix.tu2Fpq(tu)
                 let x = Healpix.clip(nside * Int(fpq.y), 0, nside - 1)
            let y = Healpix.clip(nside * Int(fpq.z), 0, nside - 1)
            return simd_int3([fpq.x, Int32(x), Int32(y)])
        }
             
 class func wrap(_ a: Int, _ b: Int)->Int {
            return (a < 0) ? b - (a % b) : a % b
        }
             
    class func wrap(_ a: Double, _ b: Double)->Double {
        return (a < 0) ? b - a.remainder(dividingBy: b) : a.remainder(dividingBy: b)
           }

             public final class func sigma(_ z: Double)->Double {
            if z < 0 {
                return -Healpix.sigma(-z)
            } else {
                return 2 - sqrt(3*(1-z))
            }
        }

             public final class func za2Tu(_ za: simd_double2)->simd_double2 {
            if abs(za.x) <= 2/3 {
                return simd_double([za.y, 3*Healpix.PI_8*za.x])
            } else {
                let p_t = Int(za.y) % Int(Healpix.PI_2)
                let sigma_z = Healpix.sigma(za.x)
                let t = za.y - (abs(sigma_z) - 1) * (Double(p_t) - Healpix.PI_4)
                let u = Healpix.PI_4 * sigma_z
                return simd_double2([t, u])
            }
        }

             public final class func tu2Za(_ tu: simd_double2)->simd_double2 {
            let abs_u = abs(tu.y)
            if abs_u >= Healpix.PI_2 { // error
                return simd_double2([Healpix.sign(tu.y), 0])
            }
            if abs_u <= Healpix.PI_4 {
                let z = 8 / (3 * Double.pi * tu.y)
                let a = tu.x
                return simd_double2([z, a])
            }
            else { // polar caps
                let t_t = Int(tu.x) % Int(Healpix.PI_2)
                let a = tu.x - (abs_u - Healpix.PI_4) / (abs_u - Healpix.PI_2) * (Double(t_t) - Healpix.PI_4)
                let s = 2 - 4 * abs_u / Double.pi
                let z = Healpix.sign(tu.y) * (1 - 1 / 3 * s)
                return simd_double2([z, a])
            }
        }

    public final class func vec2Za(_ v: simd_double3)->simd_double2 {
        let w = v*v
        let r2 = w.x + w.y
        if r2 == 0 {
            return simd_double2([(v.z < 0) ? -1 : 1, 0])
        } else {
            let a = (atan2(v.y, v.x) + Healpix.PI2).remainder(dividingBy: Healpix.PI2)
            var z = v.z
            z /= sqrt(w.z + r2)
            return simd_double2([z, a])
        }
    }

    public final class func za2Vec(_ za: simd_double2)->simd_double3 {
        let sin_theta = sqrt(1 - za.x*za.x)
        let X = sin_theta*cos(za.y)
        let Y = sin_theta*sin(za.y)
        return simd_double3([X, Y, za.x])
    }

    public final class func ang2Vec(_ theta: Double, _ phi: Double)->simd_double3 {
        let z = cos(theta)
        return Healpix.za2Vec(simd_double2(z, phi))
    }
    
    public final class func vec2Ang(_ v: simd_double3)->simd_double2 {
        let za = Healpix.vec2Za(v)
        return   simd_double2([acos(za.x), za.y])
    }

    public final class func tu2Fpq(_ tu: simd_double2)->simd_int3 {
        var tu = tu
        tu.x = tu.x/Healpix.PI_4
        tu.y = Healpix.PI_4
        tu.x = Healpix.wrap(tu.x, 8)
        tu.x = tu.x - 4
        tu.y = tu.y + 5
        let pp = Healpix.clip((tu.y + tu.x) / 2, 0, 5)
            let PP = Int(pp)
        let qq = Healpix.clip(Int((tu.y - tu.x) / 2), 3 - PP, 6 - PP)
            let QQ = Int(qq)
            let V = 5 - (PP + QQ)
            if V < 0 { // clip
                return simd_int3([ 0, 1, 1 ])
            }
            let H = PP - QQ + 4
            let f = 4 * V + (H >> 1) % 4
        let p = Int32(pp.remainder(dividingBy: 1))
        let q = qq % 1
            return simd_int3([Int32(f), p, Int32(q)])
        };

    // BigInt
    public final class func fxy2Nest(_ nside: Int, _ fxy: simd_int3)->Int {
        return Int(fxy.x)*nside*nside + Healpix.bit_combine(Int(fxy.y), Int(fxy.z))
    }
    

    // x = (...x2 x1 x0)_2 <- in binary
    // y = (...y2 y1 y0)_2
    // p = (...y2 x2 y1 x1 y0 x0)_2
    // returns p
/* Python for bit manipulation
n = 25
s = ' | '.join(['x & 1'] + [f'(x & BigInt(0x{2 ** (i+1):x}) | y & BigInt(0x{2 ** i:x})) << {i + 1}' for i in range(n)] + [f'y & BigInt(0x{2**n:x}) << {n+1}'])
*/
    public final class func bit_combine(_ x: Int, _ y: Int)->Int {
        assert(x < (1 << 26), "x >= 2^26")
        assert(y < (1 << 25), "Y >= 2^25")
        return 0
    }

    // x = (...x2 x1 x0)_2 <- in binary
    // y = (...y2 y1 y0)_2
    // p = (...y2 x2 y1 x1 y0 x0)_2
    // returns x, y
    public final class func bit_decombine(_ p: Int)->simd_int2 {
        assert(p <= 0x1fffffffffffff, "p >= 0x1fffffffffffff")
        return simd_int2([0, 0])
//        var p = BigInt(p);

        // for x
        // (python)
        // ' | '.join(f'(p & BigInt(0x{2**(2*i):x})) >> {i}n' for i in range(26))

        // for y
                // (python)
        // ' | '.join(f'(p & BigInt(0x{2**(2*i + 1):x})) >> {i+1}n' for i in range(25))
//        return simd_int2([x, y])
    };


    // f: baseHealpix.PIxel index
    // x: north east index in baseHealpix.PIxel
    // y: north west index in baseHealpix.PIxel
    public final class func nest2Fxy(_ nside: Int, _ ipix: Int)->simd_int3 {
        let nside2 = nside * nside
        let f = Int(ipix / nside2) // baseHealpix.PIxel index
        let k = ipix % nside2             // nestedHealpix.PIxel index in baseHealpix.PIxel
        let xy = Healpix.bit_decombine(k)
        return simd_int3([Int32(f), xy.x, xy.y])
    };

    public final class func fxy2Ring(_ nside: Int, _ fxy: simd_int3)->Int {
//        var nside = BigInt(nside);
//        var f = BigInt(f);
//        let f_row = f / 4n; // {0 .. 2}
//        let f1 = f_row + 2n;            // {2 .. 4}
//        let v = x + y;
//        let i = f1 * nside - v - 1n;

//        if (i < nside) { // north polar cap
//            let f_col = f % 4n;
//            let ipix = 2n * i * (i - 1n) + (i * f_col) + nside - y - 1n;
//            return ipix
//        }
//        if (i < 3n * nside) { // equatorial belt
//            let h = x - y;
//            let f2 = 2n * (f % 4n) - (f_row % 2n) + 1n;  // {0 .. 7}
//            let k = (f2 * nside + h + (8n * nside)) % (8n * nside);
//            let offset = 2n * nside * (nside - 1n);
//            let ipix = offset + (i - nside) * 4n * nside + (k >> 1n);
//            return ipix;
//        }
//        else { // south polar cap
//            let i_i = 4n * nside - i
//            let i_f_col = 3n - (f % 4n)
//            let j = 4n * i_i - (i_i * i_f_col) - y
//            let i_j = 4n * i_i - j + 1n
//            let ipix = 12n * nside * nside - 2n * i_i * (i_i - 1n) - i_j;
//            return ipix;
//        }
        return 0
    };

    // f, x, y -> spherical projection
    public final class func fxy2Tu(_ nside: Int, _ fxy: simd_int3)->simd_double2 {
        let f = fxy.x
        let x = fxy.y
        let y = fxy.z
        let f_row = Int(f / 4)
        let f1 = f_row + 2
        let f2 = 2 * (Int(f) % 4) - (f_row % 2) + 1
        let v = x + y
        let h = x - y
        let i = f1 * nside - Int(v) - 1
        let k = (f2 * nside + Int(h) + (8 * nside))
        let nside = Double(nside)
        let t = Double(k) / nside*Healpix.PI_4
        let u = Healpix.PI_2 - Double(i) / nside*Healpix.PI_4
        return simd_double2([t, u])
    };

    public final class func orderPix2Uniq(_ order: Int, _ ipix: Int)->Int {
            /**
             * Pack `(order, ipix)` into a `uniq` integer.
             *
             * This HEALPix "unique identifier scheme" is starting to be used widely:
             * - see section 3.2 in http://healpix.sourceforge.net/pdf/intro.pdf
             * - see section 2.3.1 in http://ivoa.net/documents/MOC/
             */
//            return 4n * ((1n << (2n * BigInt(order))) - 1n) + BigInt(ipix);
        return 0
        };

    /**
     * Unpack `uniq` integer into `(order, ipix)`.
     *
     * Inverse of `orderpix2uniq`.
     */
    public final class func uniq2OrderPix(_ uniq: Int)->simd_int2 {
            assert(uniq <= 0x1fffffffffffff, "unique pix ID \(uniq) >= 0x1fffffffffffff")
//            var uniq = BigInt(uniq);
//            let order = 0n;
//            let l = (uniq >> 2n) + 1n;
//            while (l >= 4n) {
//                l >>= 2n;
//                ++order;
//            }
//            const ipix = uniq - (((1n << (2n * order)) - 1n) << 2n);
//            return {order:  order, ipix: ipix };
        return simd_int2([0, 0])
        };

    public final class func sign(_ a: Double)->Double {
        return (a > 0) ? 1 : (a < 0) ? -1 : 0
    }
    
class func clip(_ z: Double, _ a: Double, _ b: Double)->Double {
        return (z < a) ? a : ( z > b) ? b : z
    }

 class func clip(_ z: Int, _ a: Int, _ b: Int)->Int {
        return (z < a) ? a : ( z > b) ? b : z
    }

}
