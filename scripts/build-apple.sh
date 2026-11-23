#!/usr/bin/env bash
# Build RatMathCore.xcframework (static library + C header) for iOS devices,
# the iOS simulator (arm64 + x86_64) and macOS. Run on a Mac with Xcode and
# rustup. Not yet run in CI; see README "Apple".
set -euo pipefail
cd "$(dirname "$0")/.."

targets=(aarch64-apple-ios aarch64-apple-ios-sim x86_64-apple-ios aarch64-apple-darwin x86_64-apple-darwin)
rustup target add "${targets[@]}"
for t in "${targets[@]}"; do cargo build --release --target "$t"; done

out=target/apple
rm -rf "$out" && mkdir -p "$out/headers" "$out/ios-sim" "$out/macos"
cp include/ratmath.h include/module.modulemap "$out/headers/"
lipo -create target/aarch64-apple-ios-sim/release/libratmath.a target/x86_64-apple-ios/release/libratmath.a -output "$out/ios-sim/libratmath.a"
lipo -create target/aarch64-apple-darwin/release/libratmath.a target/x86_64-apple-darwin/release/libratmath.a -output "$out/macos/libratmath.a"

rm -rf RatMathCore.xcframework
xcodebuild -create-xcframework \
  -library target/aarch64-apple-ios/release/libratmath.a -headers "$out/headers" \
  -library "$out/ios-sim/libratmath.a" -headers "$out/headers" \
  -library "$out/macos/libratmath.a" -headers "$out/headers" \
  -output RatMathCore.xcframework
echo "built RatMathCore.xcframework"
