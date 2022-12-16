// swift-tools-version: 5.7
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "SwiftHealpix",
    products: [
        // Products define the executables and libraries a package produces, and make them visible to other packages.
        .library(
            name: "SwiftHealpix",
            targets: ["SwiftHealpix"]),
    ],
    dependencies: [
         .package(url: "https://github.com/attaswift/BigInt.git", from: "5.3.0"),
    ],
    targets: [
        .target(
            name: "SwiftHealpix",
            dependencies: ["BigInt"]),
        .testTarget(
            name: "SwiftHealpixTests",
            dependencies: ["SwiftHealpix", "BigInt"]),
    ]
)
