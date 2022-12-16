//
//  SwiftHealpixTests.swift
//  SwiftHealpixTests
//
//  Created by Yuma decaux on 10/3/2022.
//

import XCTest
@testable import SwiftHealpix
import simd

class SwiftHealpixTests: XCTestCase {
    var hpx:Healpix!
    var tests:HealpixTests?
    
    override func setUpWithError() throws {
        try super.setUpWithError()
        hpx = Healpix(1)
        do {
            let data = try Data(contentsOf: Bundle.main.url(forResource: "tests", withExtension: "json")!)
            tests = try JSONDecoder().decode(HealpixTests.self, from: data)
        } catch let error {
            print(error.localizedDescription)
        }
    }

    override func tearDownWithError() throws {
        hpx = nil
        try super.tearDownWithError()
    }

    func testvec2Pix_nest() async {
    for test in tests!.vec2pix_nest {
        let args = test.args
        let vec = Array(args[1...args.count-1])
        let value = Healpix.vec2Pix_nest(Int(args[0])!, simd_double3(vec.map {Double($0)!}))
    let expected = Int(test.expected)!
    XCTAssertEqual(value, expected)
    }
    }

    func testvec2Pix_ring() async {
    for test in tests!.vec2pix_ring {
        let args = test.args
        let vec = Array(args[1...args.count-1])
        let value = Healpix.vec2Pix_ring(Int(args[0])!, simd_double3(vec.map {Double($0)!} ))
    let expected = Int(test.expected)!
    XCTAssertEqual(value, expected)
    }
    }

    func testang2Pix_nest() async {
    for test in tests!.ang2pix_nest {
        let args = test.args
        let value = Healpix.ang2Pix_nest(Int(args[0])!, Double(args[1])!, Double(args[2])!)
        
    let expected = Int(test.expected)!
    XCTAssertEqual(value, expected)
    }
    }

    func testang2Pix_ring() async {
    for test in tests!.ang2pix_ring {
        let args = test.args
        let value = Healpix.ang2Pix_ring(Int(args[0])!, Double(args[1])!, Double(args[2])!)
    let expected = Int(test.expected)!
    XCTAssertEqual(value, expected)
    }
    }

    func testnest2Ring() async {
    for test in tests!.nest2ring {
        let args = test.args
        let value = Healpix.nest2Ring(Int(args[0])!, Int(args[1])!)
    let expected = Int(test.expected)!
    XCTAssertEqual(value, expected)
    }
    }

    func testring2Nest() async {
    for test in tests!.ring2nest {
        let args = test.args
        let value = Healpix.ring2Nest(Int(args[0])!, Int(args[1])!)
    let expected = Int(test.expected)!
    XCTAssertEqual(value, expected)
    }
    }

    func testpix2Vec_nest() async {
    for test in tests!.pix2vec_nest {
        let args = test.args
        let value = Healpix.pix2Vec_nest(Int(args[0])!, Int(args[1])!)
        let expected = simd_double3(test.expected.components(separatedBy: ", ").map {Double($0)!} )
    XCTAssertEqual(value, expected)
    }
    }

    func testpix2Vec_ring() async {
    for test in tests!.pix2vec_ring {
        let args = test.args
        let value = Healpix.pix2Vec_ring(Int(args[0])!, Int(args[1])!)
        let expected = simd_double3(test.expected.components(separatedBy: ",").map {Double($0)!})
    XCTAssertEqual(value, expected)
    }
    }

    func testnside2PixArea() async {
    for test in tests!.nside2pixarea {
        let args = test.args
        let value = Healpix.nside2PixArea(Int(args[0])!)
    let expected = Double(test.expected)!
    XCTAssertEqual(value, expected)
    }
    }

    func testnside2Resol() async {
    for test in tests!.nside2resol {
        let args = test.args
        let value = Healpix.nside2Resol(Int(args[0])!)
    let expected = Double(test.expected)!
    XCTAssertEqual(value, expected)
    }
    }

    func testmax_pixrad() async {
    for test in tests!.max_pixrad {
        let args = test.args
        let value = Healpix.max_pixrad(Int(args[0])!)
    let expected = Double(test.expected)!
    XCTAssertEqual(value, expected)
    }
    }

    func testcorners_nest() async {
    for test in tests!.corners_nest {
        let args = test.args
        let value = Healpix.corners_nest(Int(args[0])!, Int(args[1])!)
        let expected = [simd_int3(test.expected.components(separatedBy: ", ").map {Int32($0)!})]
    XCTAssertEqual(value, expected)
    }
    }

    func testcorners_ring() async {
    for test in tests!.corners_ring {
        let args = test.args
        let value = Healpix.corners_ring(Int(args[0])!, Int(args[1])!)
        let expected = [simd_int3(test.expected.components(separatedBy: ", ").map {Int32($0)!} )]
    XCTAssertEqual(value, expected)
    }
    }

    func testorderPix2Uniq() async {
    for test in tests!.orderpix2uniq {
        let args = test.args
        let value = Healpix.orderPix2Uniq(Int(args[0])!, Int(args[1])!)
    let expected = Int(test.expected)!
    XCTAssertEqual(value, expected)
    }
    }

    func testuniq2OrderPix() async {
    for test in tests!.uniq2orderpix {
        let args = test.args
        let value = Healpix.uniq2OrderPix(Int(args[0])!)
        let expected = simd_int2(test.expected.components(separatedBy: ", ").map {Int32($0)!})
    XCTAssertEqual(value, expected)
    }
    }

}

struct HealpixTests:Decodable {
    let vec2pix_nest:[TestCase]
    let vec2pix_ring:[TestCase]
    let ang2pix_nest:[TestCase]
    let ang2pix_ring:[TestCase]
    let nest2ring:[TestCase]
    let ring2nest:[TestCase]
    let pix2vec_nest:[TestCase]
    let pix2vec_ring:[TestCase]
    let nside2pixarea:[TestCase]
    let nside2resol:[TestCase]
    let max_pixrad:[TestCase]
    let corners_nest:[TestCase]
    let corners_ring:[TestCase]
    let orderpix2uniq:[TestCase]
    let uniq2orderpix:[TestCase]
    }

struct TestCase:Decodable {
    let args:[String]
    let expected:String
}
